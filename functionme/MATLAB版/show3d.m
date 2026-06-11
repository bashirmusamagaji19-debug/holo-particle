function show3d(obj_vol, alpham, N)
%SHOW3D 旧版独立3D体绘制函数
%
%   【用途】
%   对3D体数据进行可视化渲染，基于vol3d纹理映射实现伪体绘制。
%   可调节透明度（alphamap），适用于全息重建结果的3D展示。
%
%   【在工程中的位置】
%   此为旧版独立3D显示函数。当前主GUI中嵌套了另一个show3d版本，
%   本文件作为独立调用版本保留，主要被true_way下的调试脚本使用。
%   主链路不直接调用此函数。
%
%   【输入参数】
%       obj_vol : 3D体数据矩阵（如重建强度体积）
%       alpham  : 透明度调节参数（传递给alphamap的'decrease'次数）
%       N       : 坐标轴范围尺寸（显示区间为 [0, N]）

    % 设置坐标轴刻度：0到N之间均匀取4段
    axes('fontsize', 12, 'xtick', 0:floor(N/4):N, 'ytick', 0:floor(N/4):N, 'ztick', 0:floor(N/4):N);
    % 调用vol3d进行3D纹理映射渲染，'3D'模式同时渲染三个正交方向切片
    vol3d('CData', obj_vol, 'texture', '3D');
    % 设置视角方位角-202°、仰角18°
    view(-202, 18);
    axis tight
    colormap('hot'), 
    colormap(1 - colormap);  % 反转热图色图：使高值区域显示为暗色、低值为亮色，突出粒子
    colorbar
    xlabel('z', 'fontsize', 16)
    ylabel('x', 'fontsize', 16)
    zlabel('y', 'fontsize', 16)
    axis([0 N 0 N 0 N])  % 固定坐标范围为 [0, N]
    box on,
    % zoom(0.7)
    % 降低透明度：'decrease'操作重复alpham次，使低强度区域更透明
    alphamap('decrease', alpham);
%     title('Object plot');
end
