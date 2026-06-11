function holoPre = AS_preprocessHologram(dataImg, bgImg, options)
%AS_PREPROCESSHOLOGRAM 全息图预处理（支持实验/仿真双模式）
%
%   【在整体链路中的位置】
%   本函数在 VolumeGUI_AngularSpectrum.m 的 onReconstruct 回调中
%   被调用，用于在重建之前对全息图进行必要的归一化和背景扣除。
%
%   【双模式行为】
%   - 实验模式（默认）：double 转换 → 减背景（如有）→ 输出
%   - 仿真模式（options.isSimulation = true）：
%     double 转换 → /255 归一化到 [0,1] → 减背景（如有）→ 减均值去 DC → 输出
%
%   输入：
%     dataImg  - 全息图，uint8 矩阵（来自 imread）或 double 矩阵
%     bgImg    - 背景图，可选（传 [] 表示无背景）
%     options  - 可选结构体，字段：
%                isSimulation - logical，是否为仿真全息图（默认 false）
%
%   输出：
%     holoPre  - 预处理后的全息图，double 矩阵
%                实验模式：值域约为 [0, 255]
%                仿真模式：值域约为 [-0.5, 0.5]（零均值）

    % 参数默认值
    if nargin < 3 || isempty(options)
        options = struct();
    end
    if ~isfield(options, 'isSimulation') || ~isscalar(options.isSimulation) || ~islogical(options.isSimulation)
        if isfield(options, 'isSimulation') && ~islogical(options.isSimulation)
            options.isSimulation = logical(options.isSimulation);
        else
            options.isSimulation = false;
        end
    end

    % 将输入全息图转为 double，确保后续数值运算精度
    holo = double(dataImg);

    % --- 仿真模式：先归一化到 [0, 1] ---
    % 仿真全息图以 8-bit BMP 保存 → imread 读到 uint8 [0, 255]
    % 需要恢复到 [0, 1] 区间以匹配原始的物理强度范围
    if options.isSimulation
        holo = holo / 255.0;
    end

    % --- 背景扣除（双端共用）---
    if ~isempty(bgImg)
        bg = double(bgImg);
        if options.isSimulation
            bg = bg / 255.0;
        end
        % 如果背景图与数据图尺寸不一致，强制缩放到数据图尺寸
        if ~isequal(size(holo), size(bg))
            bg = imresize(bg, size(holo));
        end
        % 直接减法：holo - bg
        holo = holo - bg;
    end

    % --- 仿真模式：去除直流分量 ---
    % 仿真全息图的 DC 被 _to_uint8 压缩到 ~0.5（归一化后）
    % 减去均值消除 FFT 零频尖峰，避免重建背景过强
    if options.isSimulation
        holo = holo - mean(holo(:));
    end

    holoPre = holo;
end
