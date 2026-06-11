function vol = AS_reconstructCS(dataImg, bgImg, ...
                                zMin_mm, zMax_mm, dz_um, ...
                                lam_nm, pix_um, ...
                                pad_size, tau, piter, tolA, iterations, useGPU)
%AS_RECONSTRUCTCS  压缩感知(CS)三维全息重建的顶层入口函数。
%
% 【在CS链路中的位置】
%   本函数是方法2（压缩感知重建）的主入口，输出完整3D幅值体。
%   与方法1（角谱快速统计，不保存完整3D体）并列。
%
% 【CS链路核心流程】
%   AS_reconstructCS（本函数）
%     → MyMakingPhase3D：生成Phase3D/Pupil（频域传播相位和光瞳）
%     → MyFieldsPropagation：生成照明场E
%     → 构造正向算子A（MyForwardOperatorPropagation → MyForwardPropagation）
%     → 构造伴随算子AT（MyAdjointOperatorPropagation → MyAdjointPropagation）
%     → 构造TV正则化 Psi/Phi
%     → 调用TwIST求解：min 0.5||y-Ax||² + tau·phi(x)
%     → MyV2C恢复复数体 → 取幅值 → 归一化得到3D体
%
% 【数据流转换】
%   TwIST只接受实向量 → MyC2V将复数展成[real;imag] → MyV2C还原
%
% 【输入】
%   dataImg    - 2D全息图矩阵（原始测量数据）
%   bgImg      - 2D背景图矩阵（用于背景扣除，可为空[]）
%   zMin_mm    - 重建起始深度 [mm]
%   zMax_mm    - 重建终止深度 [mm]
%   dz_um      - 深度层间距 [μm]
%   lam_nm     - 波长 [nm]
%   pix_um     - 探测器像素尺寸 [μm]
%   pad_size   - 零填充边宽（像素），增大频域分辨率
%   tau        - TwIST正则化参数，控制稀疏性强度
%   piter      - TV去噪内层迭代次数
%   tolA       - TwIST收敛容差
%   iterations - TwIST最大/最小迭代次数
%   useGPU     - 是否启用GPU加速（可选，默认false）
%
% 【输出】
%   vol        - 3D幅值体 [nx×ny×nz]，值域[0,1]（已归一化）
%                若有零填充，边缘已被裁剪回原始尺寸
%
% 【GPU混合模式】
%   当useGPU=true时采用hybrid策略：A/AT传播在GPU上执行（数据并行度高），
%   TV收缩和TwIST迭代在CPU上执行（避免4D TV临时缓冲区耗尽GPU显存）。

    if nargin < 13 || isempty(useGPU)
        useGPU = false;
    end
    requestedGPU = logical(useGPU);
    gpuAvailable = iHasUsableGpu();  % 检测GPU是否可用
    useGPU = requestedGPU && gpuAvailable;
    tStart = tic;

    if requestedGPU && ~gpuAvailable
        fprintf('[压缩感知] 已请求GPU，但当前不可用，改为CPU计算。\n');
    end

    try
        % 尝试用指定模式（CPU或GPU-hybrid）执行重建
        vol = iReconstructCore(dataImg, bgImg, ...
                               zMin_mm, zMax_mm, dz_um, ...
                               lam_nm, pix_um, ...
                               pad_size, tau, piter, tolA, iterations, useGPU);
        if useGPU
            modeMsg = 'GPU-hybrid';
        else
            modeMsg = 'CPU';
        end
    catch ME
        % GPU失败时自动回退到CPU模式
        if useGPU
            warning('AS_reconstructCS:GPUFallback', ...
                'GPU reconstruction failed, fallback to CPU. Reason: %s', ME.message);
            vol = iReconstructCore(dataImg, bgImg, ...
                                   zMin_mm, zMax_mm, dz_um, ...
                                   lam_nm, pix_um, ...
                                   pad_size, tau, piter, tolA, iterations, false);
            modeMsg = 'CPU (GPU fallback)';
        else
            rethrow(ME);
        end
    end

    fprintf('[压缩感知] 重建成功。模式=%s，体数据=%dx%dx%d，用时=%.2f秒\n', ...
        modeMsg, size(vol, 1), size(vol, 2), size(vol, 3), toc(tStart));
end

function vol = iReconstructCore(dataImg, bgImg, ...
                                zMin_mm, zMax_mm, dz_um, ...
                                lam_nm, pix_um, ...
                                pad_size, tau, piter, tolA, iterations, useGPU)
% iReconstructCore  重建核心逻辑（CPU或GPU-hybrid）
% ---------------------------------------------------------------
% 整体流程：
%   1. 输入验证与参数转换
%   2. 背景扣除与零填充
%   3. 生成3D传播相位与照明场
%   4. 构造正向/伴随算子、TV正则化算子
%   5. 调用TwIST求解
%   6. 后处理：复数还原 → 幅值 → 裁剪 → 归一化

    % ====== 1. 输入验证 ======
    if ~isscalar(dz_um) || ~isfinite(dz_um) || dz_um <= 0
        error('AS_reconstructCS:InvalidInput', 'dz_um must be a positive scalar.');
    end
    if ~isscalar(zMin_mm) || ~isscalar(zMax_mm) || ~isfinite(zMin_mm) || ~isfinite(zMax_mm) || zMax_mm < zMin_mm
        error('AS_reconstructCS:InvalidInput', 'zMax_mm must be greater than or equal to zMin_mm.');
    end
    if ~isscalar(lam_nm) || ~isfinite(lam_nm) || lam_nm <= 0
        error('AS_reconstructCS:InvalidInput', 'lam_nm must be a positive scalar.');
    end
    if ~isscalar(pix_um) || ~isfinite(pix_um) || pix_um <= 0
        error('AS_reconstructCS:InvalidInput', 'pix_um must be a positive scalar.');
    end
    if ~ismatrix(dataImg) || isempty(dataImg)
        error('AS_reconstructCS:InvalidInput', 'dataImg must be a non-empty 2D array.');
    end

    % 添加子函数目录到搜索路径
    if exist('./Functions', 'dir')
        addpath('./Functions');
    end

    % ====== 2. 计算精度设置 ======
    % GPU模式用single（节省显存），CPU模式用double（精度优先）
    % TwIST/TV收缩始终在CPU上运行，只有A/AT传播在GPU上
    if useGPU
        computeClass = 'single';
    else
        computeClass = 'double';
    end

    % ====== 3. 数据预处理 ======
    holo = cast(dataImg, computeClass);  % 转换计算精度
    if ~isempty(bgImg)
        bgImg = cast(bgImg, computeClass);
        if ~isequal(size(holo), size(bgImg))
            bgImg = imresize(bgImg, size(holo));  % 尺寸不匹配时自动缩放
            bgImg = cast(bgImg, computeClass);
        end
        holo = holo - bgImg;  % 背景扣除
    end

    pad_size = max(0, round(pad_size));  % 确保填充尺寸为非负整数

    % ====== 4. 物理参数转换 ======
    pixel_num = size(holo, 1);       % 原始像素数（行方向）
    detector_size = pix_um;          % 探测器像素尺寸 [μm]
    lambda_um = lam_nm / 1000;       % 波长 nm → μm
    deltaZ = dz_um;                  % 层间距 [μm]
    offsetZ = zMin_mm * 1000;        % 首层深度 mm → μm

    totalDepth_um = max(0, (zMax_mm - zMin_mm) * 1000);  % 总深度 [μm]
    nz = max(1, floor(totalDepth_um / deltaZ) + 1);       % 深度层数

    g = holo;  % 观测全息图
    shrinkage_factor = pixel_num / size(g, 1);  % 像素缩放因子（通常为1）
    sensor_size = pixel_num * detector_size;     % 传感器物理尺寸 [μm]
    deltaX = detector_size * shrinkage_factor;   % 有效x方向像素间距 [μm]
    deltaY = detector_size * shrinkage_factor;   % 有效y方向像素间距 [μm]

    % ====== 5. 零填充（增大频域分辨率） ======
    if pad_size > 0
        g = padarray(g, [pad_size pad_size]);
    end
    if useGPU
        gpuInfo = iGetGpuInfo();
        fprintf(['[压缩感知][GPU] 设备=%s，可用显存=%.2f GB，计算精度=%s，', ...
                 '模式=hybrid(A/AT on GPU, TV/TwIST on CPU)。\n'], ...
            gpuInfo.Name, gpuInfo.AvailableMemoryGB, computeClass);
    end
    [nx, ny] = size(g);  % 填充后的网格尺寸

    % ====== 6. TV正则化的3D体积尺寸参数 ======
    % MyTVpsi/MyTVphi使用(Nx,Ny,Nz)将1D实向量解释为3D数组
    % Nz=1 因为实向量维度 = 2·nx·ny·nz，被TV视为 Nx × Ny×nz·2 × 1 的3D体
    Nx = nx;
    Ny = ny * nz * 2;
    Nz = 1;

    % ====== 7. 生成照明场 ======
    % E0为均匀平面波（全1矩阵），作为MyFieldsPropagation的入射场
    if useGPU
        E0 = gpuArray.ones(nx, ny, computeClass);
    else
        E0 = ones(nx, ny, computeClass);
    end

    % 生成3D频域传播相位和光瞳函数
    [Phase3D, Pupil] = MyMakingPhase3D(nx, ny, nz, lambda_um, ...
                                       deltaX, deltaY, deltaZ, ...
                                       offsetZ, sensor_size);
    Phase3D = cast(Phase3D, computeClass);
    Pupil = cast(Pupil, computeClass);
    if useGPU
        Phase3D = gpuArray(Phase3D);  % 传输到GPU
        Pupil = gpuArray(Pupil);
    end

    % 计算照明场E：将平面波E0传播到各深度层
    E = MyFieldsPropagation(E0, nx, ny, nz, Phase3D, Pupil);

    % ====== 8. 构造观测向量（复数→实向量） ======
    g_vec = MyC2V(g(:));  % 全息图展平为[real;imag]实向量

    % ====== 9. 构造正向算子A和伴随算子AT ======
    if useGPU
        % GPU-hybrid模式：数据传入GPU执行传播，结果回CPU做TV/TwIST
        A = @(f_twist) iRunForwardHybridGpu(f_twist, E, nx, ny, nz, Phase3D, Pupil, computeClass);
        AT = @(gmeas) iRunAdjointHybridGpu(gmeas, E, nx, ny, nz, Phase3D, Pupil, computeClass);
    else
        % CPU模式：直接调用向量化接口
        A = @(f_twist) MyForwardOperatorPropagation(f_twist, E, nx, ny, nz, Phase3D, Pupil);
        AT = @(gmeas) MyAdjointOperatorPropagation(gmeas, E, nx, ny, nz, Phase3D, Pupil);
    end

    % ====== 10. 构造TV正则化算子 ======
    % Psi: TV去噪算子（近端映射），对输入做TV收缩
    % Phi: TV正则化函数值计算
    Psi = @(f, th) MyTVpsi(f, th, 0.05, piter, Nx, Ny, Nz);
    Phi = @(f) MyTVphi(f, Nx, Ny, Nz);

    % ====== 11. 调用TwIST求解 ======
    % 求解：min 0.5||g_vec - A(f)||² + tau·Phi(f)
    % Initialization=2: 初始化为AT(g_vec)，即伴随反投影
    % Monotone=1: 强制目标函数单调下降
    % StopCriterion=1: 基于目标函数相对变化收敛
    [f_vec, ~, ~, ~, ~, ~] = ...
        TwIST(g_vec, A, tau, ...
              'AT', AT, ...
              'Psi', Psi, ...
              'Phi', Phi, ...
              'Initialization', 2, ...
              'Monotone', 1, ...
              'StopCriterion', 1, ...
              'MaxIterA', iterations, ...
              'MinIterA', iterations, ...
              'ToleranceA', tolA, ...
              'Verbose', 1);

    % ====== 12. 后处理：实向量→复数体→幅值体 ======
    f_c = MyV2C(f_vec);              % 实向量还原为复数向量
    f_c = reshape(f_c, nx, ny, nz);  % 重塑为3D复数体
    mag = abs(f_c);                   % 取幅值

    % 裁剪零填充的边缘，恢复原始尺寸
    if pad_size > 0
        mag = mag(pad_size+1:end-pad_size, ...
                  pad_size+1:end-pad_size, :);
    end

    % 归一化到[0,1]
    maxVal = max(mag(:));
    if maxVal > 0
        vol = mag / maxVal;
    else
        vol = mag;
    end
end

function tf = iHasUsableGpu()
% iHasUsableGpu  检测系统是否有可用的GPU设备
    tf = false;
    try
        if exist('gpuDeviceCount', 'file') == 2 && gpuDeviceCount > 0
            gpuDevice;
            tf = true;
        end
    catch
        tf = false;
    end
end

function info = iGetGpuInfo()
% iGetGpuInfo  获取当前GPU的设备名称和可用显存
    g = gpuDevice;
    info = struct( ...
        'Name', char(g.Name), ...
        'AvailableMemoryGB', double(g.AvailableMemory) / 1024^3);
end

function y = iRunForwardHybridGpu(f_twist, E, nx, ny, nz, Phase3D, Pupil, computeClass)
% iRunForwardHybridGpu  GPU-hybrid模式下的正向算子封装
%   将TwIST的CPU实向量传入GPU → 执行正向传播 → 结果回CPU列向量
    f_twist = iToGpuArray(f_twist, computeClass);  % CPU → GPU
    y = MyForwardOperatorPropagation(f_twist, E, nx, ny, nz, Phase3D, Pupil);
    y = iGatherColumn(y);  % GPU → CPU列向量
end

function y = iRunAdjointHybridGpu(gmeas, E, nx, ny, nz, Phase3D, Pupil, computeClass)
% iRunAdjointHybridGpu  GPU-hybrid模式下的伴随算子封装
%   将TwIST的CPU实向量传入GPU → 执行伴随传播 → 结果回CPU列向量
    gmeas = iToGpuArray(gmeas, computeClass);  % CPU → GPU
    y = MyAdjointOperatorPropagation(gmeas, E, nx, ny, nz, Phase3D, Pupil);
    y = iGatherColumn(y);  % GPU → CPU列向量
end

function y = iToGpuArray(x, computeClass)
% iToGpuArray  将数据转换为指定精度并移至GPU
    y = cast(x, computeClass);
    if ~isa(y, 'gpuArray')
        y = gpuArray(y);
    end
end

function y = iGatherColumn(x)
% iGatherColumn  将GPU数组回收到CPU并展平为列向量
    if isa(x, 'gpuArray')
        x = gather(x);
    end
    y = x(:);
end
