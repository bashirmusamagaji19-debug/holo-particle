function img = AS_normalizeImage(img)
%AS_NORMALIZEIMAGE 将任意二维图像线性归一化到 [0, 1] 区间
%
%   【在整体链路中的位置】
%   本函数是方法1链路中的通用归一化工具，在以下位置被调用：
%   1. AS_reconstructAngularSpectrumFastStats.m 的 iDetectCandidates 中：
%      对 MIP 图归一化后做二值化阈值分割
%   2. AS_reconstructAngularSpectrumFastStats.m 的 iEstimateFocusDiameterFromPatch 中：
%      对焦面局部 patch 归一化后做粒径估计
%   3. AS_localizeParticles3D.m 中：
%      对完整体 MIP 图归一化后做候选检测
%   4. VolumeGUI_AngularSpectrum.m 的逆衍射逐层测试中：
%      对每帧传播结果归一化后显示
%
%   【算法逻辑】
%   线性映射：img_out = (img - min) / (max - min)
%   - 当 max > min 时：标准线性归一化到 [0, 1]
%   - 当 max == min 时：输出全零矩阵（避免除零，同时表示该图像无有效对比度）
%
%   输入：
%     img  - 任意二维图像，任意数值类型
%            常见输入：single 或 double 型的 MIP 图、焦面 patch 等
%
%   输出：
%     img  - 归一化后的二维图像，double 型，值域 [0, 1]
%            全零矩阵表示输入图像无有效对比度
%
%   【注意】
%   - 本函数保留相对亮度关系，不改变空间结构
%   - 对含负值的输入（如减背景后的场），负值会被映射到 [0, 1] 内
%   - 如果需要对含负值的场做"保留符号的归一化"（映射到 [-1, 1]），
%     应使用按最大绝对值归一化，而非本函数

    % 转为 double 确保数值精度
    img = double(img);

    % 计算全局最小值和最大值
    imgMin = min(img(:));
    imgMax = max(img(:));

    % 线性归一化到 [0, 1]
    if imgMax > imgMin
        img = (img - imgMin) / (imgMax - imgMin);
    else
        % 图像无有效对比度（常数图），输出全零
        img = zeros(size(img));
    end
end
