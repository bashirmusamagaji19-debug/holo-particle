function vol = AS_reconstructAngularSpectrum(holoPre, zVec, lambda, pixel, useGPU)
% Reference implementation only. The active method-1 baseline for Python
% alignment is Fresnel_reconstructFastStats.m.
%AS_RECONSTRUCTANGULARSPECTRUM 完整角谱�?D重建
%   vol = AS_reconstructAngularSpectrum(holoPre, zVec, lambda, pixel, useGPU)
%
%   【用途�?%   基于角谱传播理论，对预处理后的全息图在多个深度平面进行数值重建，
%   生成3D强度体数据。支持GPU加速（自动检测与回退）�?%
%   【在工程中的位置�?%   此函数是方法1（角谱法）的完整版重建实现。当前主GUI流程已不再直接调用此函数�?%   而是�?AS_reconstructAngularSpectrumFastStats（快速统计版）�?%   本函数作为备�?参考保留，主链路：GUI �?FastStats，不走此函数�?%
%   【输入参数�?%       holoPre : 预处理后的全息图（Ny × Nx 二维矩阵�?%       zVec    : 重建深度位置向量，单位：米（1 × Nz �?Nz × 1�?%       lambda  : 波长，单位：米（正标量）
%       pixel   : 像素间距，单位：米（正标量）
%       useGPU  : 是否尝试GPU加速（true/false），可选，默认false
%
%   【输出参数�?%       vol     : 重建强度体数据（Ny × Nx × Nz，single类型�?
    %% ---- 参数默认值处�?----
    if nargin < 5 || isempty(useGPU)
        useGPU = false;  % 未指定时默认使用CPU
    end

    %% ---- 输入合法性校�?----
    if isempty(holoPre) || ~ismatrix(holoPre)
        error('AS_reconstructAngularSpectrum:InvalidInput', ...
              'holoPre must be a non-empty 2D array.');
    end
    if isempty(zVec)
        error('AS_reconstructAngularSpectrum:InvalidInput', ...
              'zVec must not be empty.');
    end
    if ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
        error('AS_reconstructAngularSpectrum:InvalidInput', ...
               'lambda must be a positive scalar.');
    end
    if ~isscalar(pixel) || ~isfinite(pixel) || pixel <= 0
        error('AS_reconstructAngularSpectrum:InvalidInput', ...
               'pixel must be a positive scalar.');
    end

    %% ---- 数据类型统一为single，节省内�?----
    holo = single(holoPre);
    zVec = single(zVec(:).');  % 确保为行向量
    lambda = single(lambda);
    pixel = single(pixel);

    [ny, nx] = size(holo);
    nz = numel(zVec);

    %% ---- GPU可用性检�?----
    requestedGPU = logical(useGPU);
    gpuAvailable = iHasUsableGpu();         % 调用内部函数检测GPU
    useGPU = requestedGPU && gpuAvailable;   % 仅当请求且可用时启用

    if requestedGPU && ~gpuAvailable
        fprintf('[角谱法] 已请求GPU，但当前不可用，改为CPU计算。\n');
    end

    %% ---- 执行重建（GPU或CPU�?----
    tAll = tic;
    if useGPU
        try
            % 尝试GPU路径
            [vol, stat] = iReconstructGpuCore(holo, zVec, lambda, pixel);
            modeMsg = 'GPU';
            fprintf(['[角谱法][GPU] 上传 %.3fs，频域准�?%.3fs，FFT %.3fs，传�?IFFT %.3fs�?, ...
                     '回传 %.3fs，归一�?%.3fs，批大小 %d。\n'], ...
                stat.uploadSec, stat.gridSec, stat.fftSec, stat.propagationSec, ...
                stat.gatherSec, stat.normalizeSec, stat.batchSize);
        catch ME
            % GPU失败时自动回退CPU
            warning('AS_reconstructAngularSpectrum:GPUFallback', ...
                'GPU重建失败，回退CPU。原�? %s', ME.message);
            [vol, stat] = iReconstructCpuCore(holo, zVec, lambda, pixel);
            modeMsg = 'CPU (GPU fallback)';
            fprintf(['[角谱法][CPU回退] 频域准备 %.3fs，FFT %.3fs，传�?IFFT %.3fs�?, ...
                     '归一�?%.3fs，批大小 %d。\n'], ...
                stat.gridSec, stat.fftSec, stat.propagationSec, ...
                stat.normalizeSec, stat.batchSize);
        end
    else
        % CPU路径
        [vol, stat] = iReconstructCpuCore(holo, zVec, lambda, pixel);
        modeMsg = 'CPU';
        fprintf(['[角谱法][CPU] 频域准备 %.3fs，FFT %.3fs，传�?IFFT %.3fs�?, ...
                 '归一�?%.3fs，批大小 %d。\n'], ...
            stat.gridSec, stat.fftSec, stat.propagationSec, ...
            stat.normalizeSec, stat.batchSize);
    end

    fprintf('[角谱法] 重建完成。模�?%s，体数据=%dx%dx%d，总耗时 %.3f 秒。\n', ...
            modeMsg, ny, nx, nz, toc(tAll));
end

%% ========================================================================
%  内部函数：GPU核心重建
%  将全息图和频域网格上传到GPU，分批执行角谱传播，自动处理显存不足
%% ========================================================================
function [vol, stat] = iReconstructGpuCore(holo, zVec, lambda, pixel)
    [ny, nx] = size(holo);
    nz = numel(zVec);
    k = single(2*pi) / lambda;           % 波数 k = 2π/λ
    iUnit = complex(single(0), single(1)); % 等价�?1i，single精度

    % 计时统计结构�?    stat = struct('uploadSec', 0, 'gridSec', 0, 'fftSec', 0, ...
                  'propagationSec', 0, 'gatherSec', 0, ...
                  'normalizeSec', 0, 'batchSize', 1);

    g = gpuDevice;
    availMem = g.AvailableMemory;  % 获取当前可用显存

    %% ---- 步骤1：数据上传到GPU ----
    tUpload = tic;
    holoG = gpuArray(holo);   % 全息图上�?    zVecG = gpuArray(zVec);   % 深度向量上传
    wait(g);                   % 等待上传完成
    stat.uploadSec = toc(tUpload);

    %% ---- 步骤2：频域网格与传递函数基础构建 ----
    tGrid = tic;
    [FX, FY] = iFrequencyGrid(nx, ny, pixel, 'single'); % 生成频域坐标网格
    FX = gpuArray(FX);
    FY = gpuArray(FY);
    lambda2 = lambda * lambda;
    % 判别式：1 - λ²(fx² + fy²)，用于判断是否在传播通带�?    radicand = single(1) - lambda2 .* (FX.^2 + FY.^2);
    passband = single(radicand >= 0);       % 通带掩模：只允许传播�?    kz = sqrt(max(radicand, single(0)));    % 纵向波数分量 kz = sqrt(radicand)
    phaseBase = (iUnit * k) .* kz;          % 传播相位基础：i·k·kz（后续乘以z即为传递函数）
    wait(g);
    stat.gridSec = toc(tGrid);

    %% ---- 步骤3：全息图FFT ----
    tFft = tic;
    U0f = fft2(holoG);   % 频域表示
    wait(g);
    stat.fftSec = toc(tFft);

    %% ---- 步骤4：分批传�?----
    batchSize = iChooseBatchSizeGPU(nx, ny, nz, availMem); % 根据显存选择批大�?    stat.batchSize = batchSize;
    storeAllOnGpu = iCanStoreVolumeOnGpu(nx, ny, nz, availMem); % 判断体数据能否全部存于GPU

    currentBatch = batchSize;
    while true
        gatherInLoopSec = 0;
        try
            % 根据显存情况选择体数据存储位�?            if storeAllOnGpu
                volG = gpuArray.zeros(ny, nx, nz, 'single'); % 全部在GPU�?                vol = [];
            else
                vol = zeros(ny, nx, nz, 'single');            % 在CPU上分�?                volG = [];
            end

            tProp = tic;
            for s = 1:currentBatch:nz
                e = min(s + currentBatch - 1, nz);
                idx = s:e;

                zBatch = reshape(zVecG(idx), 1, 1, []); % 当前批次的深度向量（1×1×batch�?                H = exp(phaseBase .* zBatch);             % 角谱传递函�?H = exp(i·k·kz·z)
                H = H .* passband;                        % 应用通带掩模，滤除倏逝波

                Uz = ifft2(U0f .* H);     % 频域相乘后IFFT回到空域
                I = abs(Uz).^2;            % 取强�?
                if storeAllOnGpu
                    volG(:, :, idx) = I;                % 直接写入GPU体数�?                else
                    tg = tic;
                    vol(:, :, idx) = gather(I);         % 从GPU回传到CPU
                    gatherInLoopSec = gatherInLoopSec + toc(tg);
                end
            end
            wait(g);
            stat.propagationSec = toc(tProp);
            stat.batchSize = currentBatch;

            %% ---- 步骤5：数据回�?----
            if storeAllOnGpu
                stat.normalizeSec = 0;
                tg = tic;
                vol = gather(volG);  % 一次性从GPU回传整个体数�?                stat.gatherSec = toc(tg);
            else
                stat.gatherSec = gatherInLoopSec;  % 已在循环中逐批回传
                stat.normalizeSec = 0;
            end
            break;  % 成功完成，退出while循环

        catch ME
            % 显存不足时自动减半批大小重试
            if iIsGpuOutOfMemory(ME) && currentBatch > 1
                currentBatch = max(1, floor(currentBatch / 2));
                fprintf('[角谱法][GPU] 显存不足，批大小调整�?%d 后重试。\n', currentBatch);
                stat.batchSize = currentBatch;
                clear vol volG;
            else
                rethrow(ME);  % 非显存错误，直接抛出
            end
        end
    end

end

%% ========================================================================
%  内部函数：CPU核心重建
%  与GPU版本逻辑一致，但在CPU上执行，无显存管�?%% ========================================================================
function [vol, stat] = iReconstructCpuCore(holo, zVec, lambda, pixel)
    [ny, nx] = size(holo);
    nz = numel(zVec);
    k = single(2*pi) / lambda;           % 波数
    iUnit = complex(single(0), single(1)); % 单位虚数

    stat = struct('uploadSec', 0, 'gridSec', 0, 'fftSec', 0, ...
                  'propagationSec', 0, 'gatherSec', 0, ...
                  'normalizeSec', 0, 'batchSize', 1);

    %% ---- 频域网格与传递函数基础 ----
    tGrid = tic;
    [FX, FY] = iFrequencyGrid(nx, ny, pixel, 'single');
    lambda2 = lambda * lambda;
    radicand = single(1) - lambda2 .* (FX.^2 + FY.^2);
    passband = single(radicand >= 0);
    kz = sqrt(max(radicand, single(0)));
    phaseBase = (iUnit * k) .* kz;       % 传播相位基础
    stat.gridSec = toc(tGrid);

    %% ---- 全息图FFT ----
    tFft = tic;
    U0f = fft2(holo);
    stat.fftSec = toc(tFft);

    %% ---- 分批传播 ----
    batchSize = iChooseBatchSizeCPU(nx, ny, nz);  % CPU批大小（基于目标内存上限�?    stat.batchSize = batchSize;
    vol = zeros(ny, nx, nz, 'single');

    tProp = tic;
    for s = 1:batchSize:nz
        e = min(s + batchSize - 1, nz);
        idx = s:e;
        zBatch = reshape(zVec(idx), 1, 1, []); % 当前批次的深�?
        H = exp(phaseBase .* zBatch);  % 角谱传递函�?        H = H .* passband;             % 通带掩模
        Uz = ifft2(U0f .* H);          % 频域传播后IFFT
        vol(:, :, idx) = abs(Uz).^2;   % 取强度存入体数据
    end
    stat.propagationSec = toc(tProp);

    stat.normalizeSec = 0;
end

%% ========================================================================
%  内部函数：生成频域坐标网�?%  返回 ifftshift 后的频域坐标，与FFT输出频率排列一�?%% ========================================================================
function [FX, FY] = iFrequencyGrid(nx, ny, pixel, dtype)
    fx = ifftshift((-floor(nx/2):ceil(nx/2)-1) ./ (double(nx) * double(pixel)));
    fy = ifftshift((-floor(ny/2):ceil(ny/2)-1) ./ (double(ny) * double(pixel)));
    [FX, FY] = meshgrid(cast(fx, dtype), cast(fy, dtype));
end

%% ========================================================================
%  内部函数：根据可用显存选择GPU批大�?%  每页内存开销 = nx*ny*(8+8+8+4)字节（复数U0f + 复数H + 复数Uz + 实数I�?%  可用显存�?2%，并预留0.8GB缓冲；大图限制批大小上限�?4
%% ========================================================================
function batchSize = iChooseBatchSizeGPU(nx, ny, nz, availMem)
    perPageBytes = double(nx) * double(ny) * (8 + 8 + 8 + 4);
    usableMem = max(0, 0.52 * double(availMem) - 0.8e9);  % 52%显存减去0.8GB缓冲
    if usableMem <= 0
        batchSize = 1;
        return;
    end

    batchSize = floor(usableMem / perPageBytes);
    if nx >= 2048 || ny >= 2048
        batchSize = min(batchSize, 24);   % 大图限制批大小，避免显存碎片
    else
        batchSize = min(batchSize, 32);
    end
    batchSize = max(1, min(nz, batchSize)); % 不超过总深度数
end

%% ========================================================================
%  内部函数：选择CPU批大�?%  目标峰值内存约1.5GB
%% ========================================================================
function batchSize = iChooseBatchSizeCPU(nx, ny, nz)
    perPageBytes = double(nx) * double(ny) * (8 + 8 + 8 + 4);
    targetMem = 1.5e9;  % 目标内存上限1.5GB
    batchSize = floor(targetMem / perPageBytes);
    batchSize = max(1, min(nz, batchSize));
end

%% ========================================================================
%  内部函数：判断体数据能否全部存于GPU显存
%  条件：可用显�?> 1.4×体数据大�?+ 1GB缓冲
%% ========================================================================
function tf = iCanStoreVolumeOnGpu(nx, ny, nz, availMem)
    volBytes = double(nx) * double(ny) * double(nz) * 4;  % single类型4字节/元素
    tf = double(availMem) > (1.4 * volBytes + 1.0e9);
end

%% ========================================================================
%  内部函数：检测GPU OOM错误
%  通过错误消息关键词判断是否为显存不足
%% ========================================================================
function tf = iIsGpuOutOfMemory(ME)
    msg = lower(ME.message);
    tf = contains(msg, 'out of memory') || contains(msg, 'insufficient memory') || ...
         contains(msg, 'not enough memory') || contains(msg, 'memory');
end

%% ========================================================================
%  内部函数：检测是否有可用的GPU
%  尝试调用gpuDeviceCount，若失败则返回false
%% ========================================================================
function tf = iHasUsableGpu()
    tf = false;
    try
        if exist('gpuDeviceCount', 'file') == 2
            tf = gpuDeviceCount > 0;  % 至少有一块GPU
        end
    catch
        tf = false;
    end
end
