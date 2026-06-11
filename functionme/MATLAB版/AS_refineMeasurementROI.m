function out = AS_refineMeasurementROI(focusSlice, centerXY, searchBox, params)
%AS_REFINEMEASUREMENTROI 在最佳焦面 patch 上细化测量 ROI
%
%   【在整体链路中的位置】
%   本函数在焦层估计完成后被调用，用于将较大的搜索 ROI 收缩到
%   粒子的主能量区域，避免用过大窗口做粒径测量导致离焦拖尾混入。
%   调用位置：
%   - 方法1：AS_reconstructAngularSpectrumFastStats.m 的 iFinalizeSummaryFromMaps 中
%   - 方法2：AS_localizeParticles3D.m 中
%
%   【核心算法思想】
%   1. 从焦面 slice 中裁出搜索区域
%   2. 用边缘中位数估计背景电平
%   3. 构造非负能量图 energyMap = max(patch - bgLevel, 0)
%   4. 以能量图质心重新估计粒子中心（比 MIP 上的 WeightedCentroid 更精准）
%   5. 从质心出发，计算累计能量达到 measureRoiEnergyFraction 的最小半径
%   6. 形成更紧凑的测量窗口 measureBox
%
%   输入：
%     focusSlice - 焦面局部强度 patch，二维矩阵
%                  由 iCollectBestPatchesCpu/Gpu 返回
%                  物理含义：候选粒子在最佳焦层上的局部传播强度
%     centerXY   - 1x2 向量 [cx, cy]，候选粒子在 focusSlice 坐标系中的中心
%                  通常来自 searchBox 内的局部坐标
%     searchBox  - 1x4 向量 [x1, y1, width, height]，搜索窗口
%                  用于从 focusSlice 中裁出工作区域
%     params     - 参数结构体，可选字段：
%                  measureRoiEnergyFraction - 测量 ROI 能量占比阈值，默认 0.92
%                                             即累计 92% 能量所对应的半径
%                  roiMinRadius - 测量 ROI 最小半径（px），默认 3
%
%   输出：
%     out - 结构体，字段：
%           measureRadiusPx    - 测量 ROI 半径（px）
%           measureBox         - 1x4 [x1, y1, width, height]，全局坐标下的测量窗口
%           measureBoxLocal    - 1x4 [x1, y1, width, height]，搜索 patch 内的局部测量窗口
%           refinedCenterXY    - 1x2 [cx, cy]，细化后的全局中心坐标
%           refinedCenterLocal - 1x2 [cx, cy]，搜索 patch 内的细化局部中心
%           backgroundLevel    - 标量，估计的背景电平

    % 参数默认值填充
    if nargin < 4 || isempty(params)
        params = struct();
    end
    params = iFillDefaults(params);

    % 输入合法性检查
    if isempty(focusSlice) || ~ismatrix(focusSlice)
        error('AS_refineMeasurementROI:InvalidImage', 'focusSlice must be a non-empty 2D image.');
    end
    if numel(centerXY) ~= 2 || numel(searchBox) ~= 4
        error('AS_refineMeasurementROI:InvalidInput', 'centerXY and searchBox must be valid.');
    end

    % ---- 第1步：从 focusSlice 中裁出搜索区域 ----
    % searchBox 的坐标可能超出 focusSlice 边界，iCropBox 会自动裁剪
    [searchPatch, clippedBox] = iCropBox(focusSlice, searchBox);

    % 将 centerXY 转换到裁剪后 patch 内的局部坐标
    cx0 = centerXY(1) - clippedBox(1) + 1;
    cy0 = centerXY(2) - clippedBox(2) + 1;
    cx0 = min(size(searchPatch, 2), max(1, cx0));
    cy0 = min(size(searchPatch, 1), max(1, cy0));

    % ---- 第2步：细化中心与估计背景 ----
    % 用能量图质心替代原始中心，背景用边缘中位数估计
    [refinedLocal, bgLevel, energyMap] = iRefineCenter(searchPatch, [cx0, cy0]);

    % ---- 第3步：根据累计能量分数确定测量半径 ----
    % 从质心出发，按距离排序，找到累计能量达到 measureRoiEnergyFraction 的最小半径
    radius = iEnergyFractionRadius(energyMap, refinedLocal, params.measureRoiEnergyFraction);

    % 限制测量半径不超过 patch 尺寸的一半
    maxRadius = max(1, floor((max(1, min(size(searchPatch)) - 1)) / 2));
    radius = max(1, min(maxRadius, max(params.roiMinRadius, radius)));

    % ---- 第4步：构造测量窗口并组装输出 ----
    % 在搜索 patch 内构造居中测量窗口
    measureBoxLocal = iCenteredBox(size(searchPatch), refinedLocal(1), refinedLocal(2), radius);

    % 将局部坐标转换回全局坐标
    refinedGlobal = [clippedBox(1) + refinedLocal(1) - 1, clippedBox(2) + refinedLocal(2) - 1];
    measureBox = [clippedBox(1) + measureBoxLocal(1) - 1, clippedBox(2) + measureBoxLocal(2) - 1, ...
                  measureBoxLocal(3), measureBoxLocal(4)];

    % 组装输出结构体
    out = struct();
    out.measureRadiusPx = radius;           % 测量 ROI 半径（px）
    out.measureBox = measureBox;            % 全局坐标下的测量窗口 [x1, y1, w, h]
    out.measureBoxLocal = measureBoxLocal;  % 局部坐标下的测量窗口 [x1, y1, w, h]
    out.refinedCenterXY = refinedGlobal;    % 细化后的全局中心 [cx, cy]
    out.refinedCenterLocal = refinedLocal;  % 细化后的局部中心 [cx, cy]
    out.backgroundLevel = bgLevel;          % 估计的背景电平
end

% ---------- 内部函数：参数默认值填充 ----------
function params = iFillDefaults(params)
    % measureRoiEnergyFraction: 测量 ROI 内累计能量占比阈值
    % 0.92 表示取覆盖 92% 正能量的最小半径
    if ~isfield(params, 'measureRoiEnergyFraction') || ~isscalar(params.measureRoiEnergyFraction) || ~isfinite(params.measureRoiEnergyFraction)
        params.measureRoiEnergyFraction = 0.92;
    end

    % roiMinRadius: 测量 ROI 的最小半径（px）
    % 此处默认值为 3（比搜索 ROI 的默认 15 小很多），
    % 因为测量 ROI 是在搜索 ROI 基础上细化得到的，不需要太大的下限
    if ~isfield(params, 'roiMinRadius') || ~isscalar(params.roiMinRadius) || ~isfinite(params.roiMinRadius)
        params.roiMinRadius = 3;
    end

    % 裁剪参数范围
    params.measureRoiEnergyFraction = min(max(params.measureRoiEnergyFraction, 0.5), 0.995);
    params.roiMinRadius = max(1, round(params.roiMinRadius));
end

% ---------- 内部函数：从图像中裁出指定窗口 ----------
%   自动处理边界溢出，返回裁剪后的 patch 和实际裁剪坐标
function [patch, box] = iCropBox(img, box)
    x1 = max(1, round(box(1)));           % 左边界，不小于 1
    y1 = max(1, round(box(2)));           % 上边界，不小于 1
    x2 = min(size(img, 2), x1 + round(box(3)) - 1);  % 右边界，不超出图像
    y2 = min(size(img, 1), y1 + round(box(4)) - 1);  % 下边界，不超出图像
    box = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];  % 实际裁剪坐标
    patch = img(y1:y2, x1:x2);           % 裁剪出的 patch
end

% ---------- 内部函数：细化中心坐标 ----------
%   用边缘中位数估计背景，扣除背景后构造正能量图，
%   再用能量图质心作为细化中心
function [centerLocal, bgLevel, energyMap] = iRefineCenter(patch, fallbackLocal)
    patch = double(patch);

    % 提取四边像素用于背景估计
    if numel(patch) == 1 || size(patch, 1) < 2 || size(patch, 2) < 2
        border = patch(:).';
    else
        border = [patch(1, :), patch(end, :), patch(2:end-1, 1).', patch(2:end-1, end).'];
    end

    % 背景电平：边缘中位数（比均值更鲁棒）
    bgLevel = median(border);

    % 构造非负能量图：扣除背景后截断负值
    energyMap = patch - bgLevel;
    energyMap(energyMap < 0) = 0;

    % 计算能量图质心作为细化中心
    total = sum(energyMap(:));
    if total <= 0
        % 如果能量图为空（整个 patch 低于背景），使用回退中心
        centerLocal = fallbackLocal;
        return;
    end

    % 质心坐标 = 加权平均坐标
    [xx, yy] = meshgrid(1:size(patch, 2), 1:size(patch, 1));
    centerLocal = [sum(xx(:) .* energyMap(:)) / total, sum(yy(:) .* energyMap(:)) / total];

    % 裁剪到 patch 范围内
    centerLocal(1) = min(size(patch, 2), max(1, centerLocal(1)));
    centerLocal(2) = min(size(patch, 1), max(1, centerLocal(2)));
end

% ---------- 内部函数：计算累计能量达到目标分数的半径 ----------
%   从质心出发，按像素到质心的距离排序，
%   找到累计能量达到 targetFraction 时的距离作为半径
function radius = iEnergyFractionRadius(energyMap, centerLocal, targetFraction)
    % 计算每个像素到质心的距离
    [xx, yy] = meshgrid(1:size(energyMap, 2), 1:size(energyMap, 1));
    rr = sqrt((xx - centerLocal(1)).^2 + (yy - centerLocal(2)).^2);

    % 按距离升序排序
    [rSorted, idx] = sort(rr(:), 'ascend');
    eSorted = energyMap(idx);

    % 计算累计能量分数
    total = sum(eSorted);
    if total <= 0
        radius = 1;
        return;
    end

    cumFrac = cumsum(eSorted) / total;

    % 找到累计分数首次达到目标的位置
    pickIdx = find(cumFrac >= targetFraction, 1, 'first');
    if isempty(pickIdx)
        radius = ceil(max(rSorted));  % 如果始终未达到，取最大距离
    else
        radius = ceil(rSorted(pickIdx));  % 向上取整
    end
end

% ---------- 内部函数：构造居中正方形窗口 ----------
function box = iCenteredBox(imSize, cx, cy, radius)
    x1 = max(1, floor(cx - radius));          % 左边界
    y1 = max(1, floor(cy - radius));          % 上边界
    x2 = min(imSize(2), ceil(cx + radius));   % 右边界
    y2 = min(imSize(1), ceil(cy + radius));   % 下边界
    box = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];
end
