function hParticleFig = AS_showParticleSpheres(coords3D, diams_um, pix_um, parentAx)
%AS_SHOWPARTICLESPHERES 稀疏粒子球体显示（true_way版本）
%
%   【用途】
%   将定位到的粒子3D坐标以放大球体的形式在3D空间中可视化渲染。
%   每个粒子以半透明球体表示，颜色自动分配，球体半径取粒子直径的放大值。
%   x/y坐标从像素转换为微米，z轴保持微米单位。
%
%   【在工程中的位置】
%   此为true_way目录下的历史版本。根目录下曾有同名文件但已不存在。
%   当前主GUI可能使用其他版本。本版本作为原型/调试版本保留。
%   主链路：GUI →（可能调用其他版本的showParticleSpheres），不一定走此文件。
%
%   【输入参数】
%       coords3D  : 粒子3D坐标矩阵（N×3），列顺序为 [x_pixel, y_pixel, z_um]
%       diams_um  : 粒子直径向量（N×1），单位：微米；为空时使用默认半径
%       pix_um    : 像素尺寸，单位：微米/像素，用于将x/y从像素转为微米
%       parentAx  : 父坐标轴句柄（可选）；为空时创建新图窗
%
%   【输出参数】
%       hParticleFig : 粒子显示图窗的句柄

    % 无粒子时直接返回空
    if isempty(coords3D)
        hParticleFig = [];
        return;
    end
    if nargin < 4
        parentAx = [];  % 未指定父坐标轴
    end

    %% ---- 坐标转换：像素 → 微米 ----
    % x和y列乘以像素尺寸，z列已经是微米单位无需转换
    coordsUm = [coords3D(:, 1) * pix_um, coords3D(:, 2) * pix_um, coords3D(:, 3)];

    %% ---- 计算球体基准半径 ----
    if isempty(diams_um)
        baseRadius = max(40, 6 * pix_um);  % 无直径信息时，取40μm或6像素对应微米的较大值
    else
        baseRadius = max([40, 6 * pix_um, median(diams_um) * 1.5]);  % 取三者的最大值
    end

    %% ---- 创建或复用图窗/坐标轴 ----
    if isempty(parentAx) || ~isgraphics(parentAx, 'axes')
        % 创建新的深色主题图窗
        hParticleFig = figure( ...
            'Name', 'Sparse Particle Field', ...
            'NumberTitle', 'off', ...
            'Color', [0.05 0.07 0.12], ...      % 深蓝黑色背景
            'Position', [1150 120 820 640]);
        ax = axes('Parent', hParticleFig);
    else
        % 复用已有坐标轴
        ax = parentAx;
        hParticleFig = ancestor(ax, 'figure');
        cla(ax);  % 清空坐标轴内容
    end
    hold(ax, 'on');

    %% ---- 绘制粒子球体 ----
    [sx, sy, sz] = sphere(20);  % 生成单位球面（20×20分辨率）
    colors = lines(max(size(coordsUm, 1), 7));  % 自动分配颜色，至少7种
    for n = 1:size(coordsUm, 1)
        % 每个粒子的显示半径：取基准半径与(1.2倍直径/2)的较大值
        radiusUm = max(baseRadius, diams_um(min(n, numel(diams_um))) * 1.2 / 2);
        surf(ax, ...
            sx * radiusUm + coordsUm(n, 1), ...  % x方向缩放+平移
            sy * radiusUm + coordsUm(n, 2), ...  % y方向缩放+平移
            sz * radiusUm + coordsUm(n, 3), ...  % z方向缩放+平移
            'FaceColor', colors(mod(n-1, size(colors, 1)) + 1, :), ... % 循环取色
            'EdgeColor', 'none', ...              % 无边框线
            'FaceAlpha', 0.95);                   % 近不透明
    end

    %% ---- 设置坐标轴外观 ----
    grid(ax, 'on');
    axis(ax, 'equal');            % 等比例显示
    view(ax, 32, 24);            % 设定视角方位角32°、仰角24°
    xlabel(ax, ['x (' char(181) 'm)']);  % μm符号
    ylabel(ax, ['y (' char(181) 'm)']);
    zlabel(ax, ['z (' char(181) 'm)']);
    title(ax, 'Localized Particle Field');
    % 深色主题配色
    set(ax, 'Color', [0.10 0.11 0.15], ...        % 坐标区背景色
        'XColor', [0.85 0.90 1.00], ...            % x轴颜色（浅蓝白）
        'YColor', [0.85 0.90 1.00], ...
        'ZColor', [0.85 0.90 1.00], ...
        'FontName', 'Times New Roman');
    camlight(ax, 'headlight');  % 头顶光源
    lighting(ax, 'gouraud');    % Gouraud平滑光照
    drawnow;
end
