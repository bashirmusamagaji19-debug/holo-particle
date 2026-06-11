function VolumeGUI_AngularSpectrum
% VolumeGUI_AngularSpectrum - 3D数字全息重建GUI主函数
%   该函数创建一个图形用户界面，用于加载全息图和背景图，
%   设置重建参数，执行角谱重建或压缩感知重建，并进行粒子统计分析。
%
%   GUI布局：
%   - 左列：重建参数设置
%   - 中列：输入全息图与背景图显示
%   - 右列：统计结果与可视化

    % 添加函数库路径
    addpath('./function/');

    %==================== 创建主窗口 ====================%
    hFig = figure( ...
        'Name', '3D Digital Hologram Reconstruction', ... % 窗口标题
        'NumberTitle', 'off', ... % 不显示默认的Figure编号
        'Color', [0.05 0.07 0.12], ... % 背景色：深蓝灰色
        'Position', [100 100 1400 600]); % 窗口位置和大小 [左, 下, 宽, 高]
    try
        set(hFig, 'Renderer', 'opengl'); % 尝试使用OpenGL渲染器以加速3D显示
    catch
    end

    % 顶部白色标题文本
    uicontrol('Parent', hFig, 'Style', 'text', ...
        'String', '3D Digital Hologram Reconstruction', ...
        'Units', 'normalized', ... % 使用归一化单位，随窗口缩放
        'Position', [0.02 0.92 0.6 0.06], ... % 位置 [左, 下, 宽, 高]
        'ForegroundColor', [0.90 0.95 1.00], ... % 前景色（文字）：亮白
        'BackgroundColor', [0.05 0.07 0.12], ... % 背景色：与窗口一致
        'FontSize', 18, ... % 字体大小
        'FontWeight', 'bold', ... % 加粗
        'HorizontalAlignment', 'left'); % 左对齐

    % 定义界面配色方案
    primaryColor   = [0.00 0.65 0.70]; % 主色调：青色
    secondaryColor = [0.16 0.32 0.63]; % 次要色：深蓝色
    accentColor    = [0.30 0.50 0.80]; % 强调色：浅蓝色
    textColor      = [0.85 0.90 1.00]; % 文字颜色：亮白
    panelBg        = [0.09 0.10 0.14]; % 面板背景色：深灰

    % 定义特殊字符（用于单位显示）
    muChar  = char(181);   % 微米符号 'µ'
    lamChar = char(955);   % 波长符号 'λ'

    % 初始化handles结构体，用于存储所有数据和控件句柄
    handles = struct();
    handles.DataImage     = []; % 存储加载的数据全息图
    handles.BgImage       = []; % 存储加载的背景图
    handles.VolumeData    = []; % 存储3D重建体数据
    handles.PreprocessedHologram = []; % 存储预处理后的全息图
    handles.VolumeZVecUm  = []; % 存储z轴坐标向量（单位：微米）
    handles.FastStatsSummary = []; % 存储快速统计模式的结果
    handles.ParticleCoords3D = []; % 存储粒子的3D坐标
    handles.ParticleAxialCurves = {}; % 存储粒子的轴向曲线
    handles.StatsResult   = []; % 存储最终的统计结果
    handles.VolumeAlpha   = 0.1;   % 3D渲染的透明度参数
    handles.ZoomApplied   = false; % 标记是否已应用相机缩放
    handles.RenderTimer   = []; % 用于异步渲染的定时器
    handles.RenderRequestId = 0; % 渲染请求ID，用于取消旧请求
    handles.RenderTitle   = ''; % 渲染标题
    handles.FastRenderEnabled = true; % 是否启用快速渲染（降低分辨率以加速）
    handles.FastRenderMaxDim = 256; % 快速渲染的最大维度
    handles.FastRenderTexture = '2D'; % 快速渲染优先使用2D纹理

    % --- [新增] 顶部菜单栏 ---
    hMenu = uimenu(hFig, 'Label', '分析 (Analysis)');
    hStatsMenu = uimenu(hMenu, ...
        'Label', '生成粒径统计报告 (Generate Particle Report)', ...
        'Callback', @onSimpleStatsV2); % 点击后调用统计函数
    try
        set(hStatsMenu, 'Interruptible', 'on', 'BusyAction', 'queue');
    catch
        set(hStatsMenu, 'Interruptible', 'on');
    end
    % ------------------------
    


    %==================== 左列：参数面板 ====================%
    leftPanel = uipanel('Parent', hFig, ...
        'Title', ' Reconstruction Parameters ', ... % 面板标题
        'Units', 'normalized', ...
        'Position', [0.01 0.10 0.18 0.8], ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', textColor, ...
        'BorderType', 'line', ...
        'HighlightColor', [0.4 0.4 0.4], ...
        'ShadowColor', [0 0 0], ...
        'FontSize', 12, 'FontWeight', 'bold');

    % --- z_min (mm) ---
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'z最小值 (mm):', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.83 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 12, ...
        'BackgroundColor', panelBg, ...
        'HorizontalAlignment', 'left');
    handles.hZmin = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.83 0.45 0.06], ...
        'String', '20'); % 默认值

    % --- z_max (mm) ---
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'z最大值 (mm):', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.75 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 12, ...
        'HorizontalAlignment', 'left');
    handles.hZmax = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.75 0.45 0.06], ...
        'String', '35'); % 默认值

    % --- dz (µm) ---
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', ['z步长 (' muChar 'm):'], ...
        'Units', 'normalized', ...
        'Position', [0.05 0.67 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hDz = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.67 0.45 0.06], ...
        'String', '50'); % 默认值

    % --- 波长 (nm) ---
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', ['波长 (nm):'], ...
        'Units', 'normalized', ...
        'Position', [0.05 0.59 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hLam = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.59 0.45 0.06], ...
        'String', '638'); % 默认值：638nm (红光)

    % --- 像素大小 (µm) -> detector_size ---
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', ['像素尺寸 (' muChar 'm):'], ...
        'Units', 'normalized', ...
        'Position', [0.05 0.51 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hPixel = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.51 0.45 0.06], ...
        'String', '3.45'); % 默认值

    % --- 重建方法选择 ---
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', '重建方法:', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.43 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hMethodPopup = uicontrol('Parent', leftPanel, 'Style', 'popupmenu', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.43 0.45 0.06], ...
        'BackgroundColor', [0.15 0.17 0.22], ...
        'ForegroundColor', textColor, ...
        'String', {'方法1：反向衍射'});

    
    % --- GPU加速复选框 ---
    handles.hUseGPU = uicontrol('Parent', leftPanel, 'Style', 'checkbox', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.01 0.88 0.04], ...
        'String', 'Use GPU acceleration (if available)', ...
        'Value', 1, ... % 默认勾选
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9);

    % --- 仿真模式复选框 ---
    handles.hSimulationMode = uicontrol('Parent', leftPanel, 'Style', 'checkbox', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.06 0.88 0.04], ...
        'String', 'Simulation mode', ...
        'Value', 0, ... % 默认不勾选（实验模式）
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9);

    %==================== 中列：数据与背景图面板 ====================%
    set(handles.hDz, 'String', '20'); % 重新设置dz默认值（可能是调试残留）

    midPanel = uipanel('Parent', hFig, ...
        'Title', ' Input Hologram & Background ', ...
        'Units', 'normalized', ...
        'Position', [0.21 0.10 0.26 0.8], ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'ForegroundColor', textColor, ...
        'BorderType', 'line', ...
        'HighlightColor', [0.4 0.4 0.4], ...
        'ShadowColor', [0 0 0], ...
        'FontSize', 12, 'FontWeight', 'bold');

    % 数据图显示轴
    handles.hDataAxes = axes('Parent', midPanel, ...
        'Units', 'normalized', ...
        'Position', [0.10 0.60 0.85 0.34]);
    title(handles.hDataAxes, '数据图');
    set(handles.hDataAxes, 'Color', [0.10 0.11 0.15], ...
        'XColor', textColor, 'YColor', textColor, 'FontName', 'Times New Roman');

    % 加载数据图按钮
    hLoadData = uicontrol('Parent', midPanel, 'Style', 'pushbutton', ...
        'String', '加载数据图', ...
        'Units', 'normalized', ...
        'Position', [0.25 0.48 0.50 0.06], ...
        'BackgroundColor', secondaryColor, ...
        'ForegroundColor', [0.95 0.97 1.00], ...
        'FontWeight', 'bold');

    % 背景图显示轴
    handles.hBgAxes = axes('Parent', midPanel, ...
        'Units', 'normalized', ...
        'Position', [0.10 0.12 0.85 0.34]);
    title(handles.hBgAxes, '');
    set(handles.hBgAxes, 'Color', [0.10 0.11 0.15], ...
        'XColor', textColor, 'YColor', textColor, 'FontName', 'Times New Roman');

    % 加载背景图按钮
    hLoadBg = uicontrol('Parent', midPanel, 'Style', 'pushbutton', ...
        'String', '加载背景图', ...
        'Units', 'normalized', ...
        'Position', [0.25 0.01 0.50 0.06], ...
        'BackgroundColor', secondaryColor, ...
        'ForegroundColor', [0.95 0.97 1.00], ...
        'FontWeight', 'bold');

    %==================== 右列：统计结果面板 ====================%
    rightPanel = uipanel('Parent', hFig, ...
        'Title', ' Particle Statistics ', ...
        'Units', 'normalized', ...
        'Position', [0.50 0.10 0.49 0.8], ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'ForegroundColor', textColor, ...
        'BorderType', 'line', ...
        'HighlightColor', [0.4 0.4 0.4], ...
        'ShadowColor', [0 0 0], ...
        'FontSize', 12, 'FontWeight', 'bold');

    % 开始计算按钮
    hReconstructPrep = uicontrol('Parent', rightPanel, 'Style', 'pushbutton', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.88 0.22 0.07], ...
        'String', '开始计算', ...
        'BackgroundColor', secondaryColor, ...
        'ForegroundColor', [0.95 0.97 1.00], ...
        'FontWeight', 'bold');

    % 导出结果按钮
    handles.hExportBtn = uicontrol('Parent', rightPanel, 'Style', 'pushbutton', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.79 0.22 0.07], ...
        'String', '导出结果', ...
        'BackgroundColor', accentColor, ...
        'ForegroundColor', [0.95 0.97 1.00], ...
        'Enable', 'off', ... % 初始禁用，有结果后启用
        'FontWeight', 'bold');
    try
        set(hReconstructPrep, 'Interruptible', 'on', 'BusyAction', 'queue');
    catch
        set(hReconstructPrep, 'Interruptible', 'on');
    end
    try
        set(handles.hExportBtn, 'Interruptible', 'on', 'BusyAction', 'queue');
    catch
        set(handles.hExportBtn, 'Interruptible', 'on');
    end
    
    % 逆衍射逐层测试按钮
    hInverseTest = uicontrol('Parent', rightPanel, 'Style', 'pushbutton', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.70 0.22 0.06], ...
        'String', '逆衍射逐层测试', ...
        'BackgroundColor', [0.92 0.78 0.34], ... % 黄色
        'ForegroundColor', [0.08 0.08 0.08], ...
        'FontWeight', 'bold');
    try
        set(hInverseTest, 'Interruptible', 'on', 'BusyAction', 'queue');
    catch
        set(hInverseTest, 'Interruptible', 'on');
    end
    
    % 状态文本显示
    handles.hStatusText = uicontrol('Parent', rightPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.58 0.30 0.10], ...
        'String', '状态：空闲', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');

    % 统计摘要面板
    summaryPanel = uipanel('Parent', rightPanel, ...
        'Title', ' 统计摘要 ', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.32 0.33 0.23], ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'ForegroundColor', textColor, ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
    handles.hSummaryText = uicontrol('Parent', summaryPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.05 0.90 0.90], ...
        'String', '尚未生成统计结果', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');

    % 粒径分布直方图轴
    handles.hHistAxes = axes('Parent', rightPanel, ...
        'Units', 'normalized', ...
        'Position', [0.40 0.56 0.56 0.34], ...
        'Color', [0.10 0.11 0.15], ...
        'XColor', textColor, 'YColor', textColor, ...
        'FontName', 'Times New Roman');
    title(handles.hHistAxes, '粒径分布', 'Color', textColor);
    xlabel(handles.hHistAxes, ['直径 (' muChar 'm)'], 'Color', textColor);
    ylabel(handles.hHistAxes, '粒子数量', 'Color', textColor);
    grid(handles.hHistAxes, 'on');

    % 粒子数据结果表格
    handles.hResultTable = uitable('Parent', rightPanel, ...
        'Units', 'normalized', ...
        'Position', [0.37 0.05 0.60 0.42], ...
        'Data', cell(0, 5), ...
        'ColumnName', {'ID', ['x (' muChar 'm)'], ['y (' muChar 'm)'], ['z (' muChar 'm)'], ['d (' muChar 'm)']}, ...
        'RowName', [], ...
        'BackgroundColor', [1 1 1], ...
        'ForegroundColor', [0 0 0]);

    % 预处理与定位参数面板
    paramPanel = uipanel('Parent', rightPanel, ...
        'Title', ' 预处理与定位参数 ', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.01 0.33 0.29], ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'ForegroundColor', textColor, ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    % 峰比值阈值
    uicontrol('Parent', paramPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.84 0.52 0.09], ...
        'String', '峰比值阈值', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');
    handles.hPeakRatio = uicontrol('Parent', paramPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.60 0.84 0.28 0.09], ...
        'String', '1.15');

    % Z边界层数
    uicontrol('Parent', paramPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.63 0.52 0.09], ...
        'String', 'Z边界层数:', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');
    handles.hEdgeMargin = uicontrol('Parent', paramPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.60 0.63 0.28 0.09], ...
        'String', '2');

    % 二维识别阈值
    uicontrol('Parent', paramPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.42 0.52 0.09], ...
        'String', '二维识别阈值', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');
    handles.hXYThreshold = uicontrol('Parent', paramPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.60 0.42 0.28 0.09], ...
        'String', '0.65');

    % ROI半径倍率
    uicontrol('Parent', paramPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.21 0.52 0.09], ...
        'String', 'ROI半径倍率:', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');
    handles.hROIScale = uicontrol('Parent', paramPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.60 0.21 0.28 0.09], ...
        'String', '0.50');

    % ROI最小半径
    uicontrol('Parent', paramPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.00 0.52 0.09], ...
        'String', 'ROI最小半径', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');
    handles.hROIMinRadius = uicontrol('Parent', paramPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.60 0.00 0.28 0.09], ...
        'String', '15');

    % 初始化快速渲染参数（如果为空）
    if ~isfield(handles, 'FastRenderMaxDim') || isempty(handles.FastRenderMaxDim)
        handles.FastRenderMaxDim = 256;
    end
    guidata(hFig, handles); % 保存handles数据

    %==================== 绑定回调函数 ====================%
    set(hLoadData,    'Callback', @onLoadData);       % 加载数据图
    set(hLoadBg,      'Callback', @onLoadBg);         % 加载背景图
    set(hReconstructPrep, 'Callback', @onReconstruct);% 开始重建
    set(handles.hExportBtn, 'Callback', @onExportResults); % 导出结果
    set(hInverseTest, 'Callback', @onPlayInverseTest);% 逆衍射测试
    set(handles.hMethodPopup, 'Callback', @onMethodChanged); % 方法切换
    handles = iRefreshUiState(handles); % 刷新UI状态（如控件启用/禁用）
    iUpdateStatsDisplay(handles); % 更新统计显示
    guidata(hFig, handles);

    %==================== 以下为嵌套回调函数 ====================%

    % --- 加载数据图回调 ---
    function onLoadData(~, ~)
        handles = guidata(hFig);
        % 打开文件选择对话框
        [file, path] = uigetfile( ...
            {'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Image Files'}, ...
            '选择数据图片(图片1)');
        if isequal(file, 0), return; end % 用户取消

        img = imread(fullfile(path, file));
        
        % 如果是RGB图，转为灰度图
        if ndims(img) == 3, img = rgb2gray(img); end
        img = double(img); % 转为双精度

        % 更新handles数据，清空旧的重建结果
        handles.DataImage  = img;
        handles.VolumeData = [];
        handles.PreprocessedHologram = [];
        handles.VolumeZVecUm = [];
        handles.FastStatsSummary = [];
        handles.ParticleCoords3D = [];
        handles.ParticleAxialCurves = {};
        handles.StatsResult = [];
        handles.ZoomApplied = false;  % 重置缩放状态
        
        axes(handles.hDataAxes);
        cla(handles.hDataAxes); % 清空坐标轴

        % 显示图像
        imagesc(handles.hDataAxes, img);
        axis(handles.hDataAxes, 'image'); % 保持纵横比
        colormap(handles.hDataAxes, 'gray'); % 灰度 colormap

        % 添加颜色条并设置字体
        cbar = colorbar(handles.hDataAxes);
        cbar.FontName = 'Times New Roman';
        cbar.FontSize = 9;
        cbar.Color    = [1 1 1];

        set(handles.hDataAxes, ...
            'FontName', 'Times New Roman', ...
            'FontSize', 9, ...
            'XColor', [1 1 1], ...
            'YColor', [1 1 1]);

        % 更新状态文本
        set(handles.hStatusText, 'String', '状态：全息图已加载');
        drawnow;

        handles = iRefreshUiState(handles);
        iUpdateStatsDisplay(handles);
        guidata(hFig, handles);
    end

    % --- 加载背景图回调 ---
    function onLoadBg(~, ~)
        handles = guidata(hFig);
        [file, path] = uigetfile( ...
            {'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Image Files'}, ...
            '选择背景图片(图片2)');
        if isequal(file, 0), return; end

        img = imread(fullfile(path, file));
        if ndims(img) == 3, img = rgb2gray(img); end
        img = double(img);

        handles.BgImage = img;
        handles.PreprocessedHologram = [];
        handles.FastStatsSummary = [];
        handles.ParticleCoords3D = [];
        handles.ParticleAxialCurves = {};
        handles.StatsResult = [];

        axes(handles.hBgAxes);
        cla(handles.hBgAxes);
        imagesc(handles.hBgAxes, img);
        axis(handles.hBgAxes, 'image');
        colormap(handles.hBgAxes, 'gray');
        colorbar(handles.hBgAxes);
        title(handles.hBgAxes, '背景图');

        set(handles.hStatusText, 'String', '状态：背景图已加载');
        drawnow;

        handles = iRefreshUiState(handles);
        iUpdateStatsDisplay(handles);
        guidata(hFig, handles);
    end

    % --- 透明度改变回调（当前界面未直接绑定此控件） ---
    function onAlphaChanged(src, ~)
        handles = guidata(hFig);
        val = get(src, 'Value');  % 获取滑块值 (0~1)
        handles.VolumeAlpha = val;
        set(handles.hAlphaText, 'String', sprintf('Alpha = %.2f', val));
        guidata(hFig, handles);
    end

    % --- 开始重建主回调 ---
    function onReconstruct(~, ~)
        handles = guidata(hFig);

        set(handles.hStatusText, 'String', '状态：正在计算粒子统计结果...');
        drawnow;

        % 检查是否加载了数据图
        if isempty(handles.DataImage)
            errordlg('请先加载数据图', '提示');
            set(handles.hStatusText, 'String', '状态：重建失败（未加载数据图像）');
            drawnow;
            return;
        end

        dataImg = handles.DataImage;

        % 处理背景图：如果有，调整大小与数据图一致
        if ~isempty(handles.BgImage)
            bgImg = handles.BgImage;
            if ~isequal(size(dataImg), size(bgImg))
                bgImg = imresize(bgImg, size(dataImg));
            end
        else
            bgImg = [];
        end
        % 读取基础重建参数
        zMin_mm = str2double(get(handles.hZmin,  'String'));
        zMax_mm = str2double(get(handles.hZmax,  'String'));
        dz_um   = str2double(get(handles.hDz,    'String'));
        lam_nm  = str2double(get(handles.hLam,   'String'));
        pix_um  = str2double(get(handles.hPixel, 'String'));

        % 检查参数有效性
        if any(isnan([zMin_mm, zMax_mm, dz_um, lam_nm, pix_um]))
            errordlg('请检查所有参数是否为数字', '参数错误');
            set(handles.hStatusText, 'String', '状态：重建失败（参数错误）');
            drawnow;
            return;
        end
        % 单位换算：转为国际单位制 (米)
        zMin   = zMin_mm * 1e-3;   % mm -> m
        zMax   = zMax_mm * 1e-3;   % mm -> m
        dz     = dz_um   * 1e-6;   % µm -> m
        lambda = lam_nm  * 1e-9;   % nm -> m
        pixel  = pix_um  * 1e-6;   % µm -> m

        % 检查z范围逻辑
        if dz <= 0 || zMax <= zMin
            errordlg('请检查 z_min、z_max 和 dz：z_max > z_min 且 dz > 0', '参数错误');
            set(handles.hStatusText, 'String', '状态：重建失败（参数无效）');
            drawnow;
            return;
        end

        % 生成z向量
        zVec = zMin:dz:zMax;
        if numel(zVec) < 2
            errordlg('z 范围过小或步长过大，导致层数不足', '参数错误');
            set(handles.hStatusText, 'String', '状态：重建失败（层数不足）');
            drawnow;
            return;
        end

        % 检查是否使用GPU
        if isfield(handles, 'hUseGPU') && isgraphics(handles.hUseGPU)
            useGPU = logical(get(handles.hUseGPU, 'Value'));
        else
            useGPU = true;
        end

        % 检查是否为仿真模式
        if isfield(handles, 'hSimulationMode') && isgraphics(handles.hSimulationMode)
            simulationMode = logical(get(handles.hSimulationMode, 'Value'));
        else
            simulationMode = false;
        end

        % --- 1. 预处理全息图 ---
        tPre = tic;
        prepOpts = struct('isSimulation', simulationMode);
        holoPre = AS_preprocessHologram(dataImg, bgImg, prepOpts);
        preSec = toc(tPre);
        if simulationMode
            fprintf('[流程耗时] 预处理完成（仿真模式：/255 归一化 + 减DC），用时 %.2f 秒。\n', preSec);
        else
            fprintf('[流程耗时] 预处理完成，用时 %.2f 秒。\n', preSec);
        end
        handles.PreprocessedHologram = holoPre;
        handles.FastStatsSummary = [];
        handles.StatsResult = [];

        % --- 2. 执行重建 ---
        tRecon = tic;
        fastStats = [];
        methodName = '方法1：反向衍射';
        % 读取粒子定位参数
        peakRatioTh = str2double(get(handles.hPeakRatio, 'String'));
        edgeMargin = str2double(get(handles.hEdgeMargin, 'String'));
        xyThreshold = str2double(get(handles.hXYThreshold, 'String'));
        roiScale = str2double(get(handles.hROIScale, 'String'));
        roiMinRadius = str2double(get(handles.hROIMinRadius, 'String'));
        if isnan(peakRatioTh), peakRatioTh = 1.15; end
        if isnan(edgeMargin), edgeMargin = 2; end
        if isnan(xyThreshold), xyThreshold = 0.85; end
        if isnan(roiScale), roiScale = 0.50; end
        if isnan(roiMinRadius), roiMinRadius = 15; end

        % 封装快速统计参数
        fastStatsParams = struct( ...
            'minPeakRatio', peakRatioTh, ...
            'edgeMargin', edgeMargin, ...
            'xyThreshold', xyThreshold, ...
            'roiScale', roiScale, ...
            'roiMinRadius', roiMinRadius, ...
            'simulationMode', simulationMode, ...
            'measureRoiEnergyFraction', 0.92, ...
            'searchRoiGrowFactor', 1.25, ...
            'edgeEnergyThreshold', 0.20, ...
            'searchMinRadiusUm', 30);

        % 调用外部函数：角谱重建+快速统计（不生成完整3D体以节省内存）
        fastStats = Fresnel_reconstructFastStats( ...
            holoPre, zVec, lambda, pixel, useGPU, fastStatsParams);
        reconSec = toc(tRecon);
        fprintf('[流程耗时] 重建完成，用时 %.2f 秒。\n', reconSec);

        % --- 3. 结果后处理 ---
        % 保存快速统计结果
        handles.VolumeData = [];
        handles.FastStatsSummary = fastStats;
        handles.VolumeZVecUm = fastStats.zVecUm(:);
        handles.ParticleCoords3D = fastStats.coords3D;
        handles.ParticleAxialCurves = fastStats.axialCurves;
        handles = iRefreshUiState(handles);
        guidata(hFig, handles);

        % --- 仿真模式：显示双阈值诊断图 ---
        if simulationMode && isfield(fastStats, 'bwLow')
            hDiag = figure('Name', '候选检测诊断 — 双阈值分割', ...
                'NumberTitle', 'off', 'Position', [150 120 1100 520]);
            % 归一化 MIP
            ax1 = subplot(2,3,1); imshow(fastStats.img2D, []); title('归一化 MIP'); colorbar;
            % 低阈值二值化
            ax2 = subplot(2,3,2); imshow(fastStats.bwLow);  title(sprintf('低阈值 (联动区域)'));
            % 高阈值二值化
            ax3 = subplot(2,3,3); imshow(fastStats.bwHigh); title(sprintf('高阈值 (亮核)'));
            % 最终 mask + 候选框
            ax4 = subplot(2,3,4);
            imshow(fastStats.img2D, []); hold on;
            visboundaries(fastStats.bw, 'Color', 'g', 'LineWidth', 0.5);
            title(sprintf('最终 Mask (%d 候选)', size(fastStats.candidateRoiBoxes, 1)));
            % 高阈值叠加在 MIP 上
            ax5 = subplot(2,3,5);
            imshow(fastStats.img2D, []); hold on;
            bwOverlay = imoverlay(mat2gray(fastStats.img2D), fastStats.bwHigh, [1 0 0]);
            imshow(bwOverlay); title('MIP + 高阈值(红)');
            % 候选框画在 MIP 上
            ax6 = subplot(2,3,6);
            imshow(fastStats.img2D, []); hold on;
            for nb = 1:size(fastStats.candidateRoiBoxes, 1)
                rectangle('Position', fastStats.candidateRoiBoxes(nb,:), ...
                    'EdgeColor', 'g', 'LineWidth', 1);
            end
            title('MIP + 候选 ROI 框');
            linkaxes([ax1 ax2 ax3 ax4 ax5 ax6]);
        end

        onSimpleStatsV2([], []); % 直接调用统计显示
        return;
    end

    % --- 逆衍射逐层测试回调 ---
    function onPlayInverseTest(~, ~)
        handles = guidata(hFig);
        set(handles.hStatusText, 'String', '状态：正在准备逆衍射逐层测试...');
        drawnow;

        % 收集测试所需的输入参数
        [ctx, errMsg] = iCollectInversePlaybackInputs(handles);
        if ~isempty(errMsg)
            errordlg(errMsg, '参数错误');
            set(handles.hStatusText, 'String', '状态：逆衍射逐层测试失败');
            drawnow;
            return;
        end

        % 检查是否为方法1（仅方法1支持此测试）
        if ~(isfield(ctx, 'methodIndex') && ctx.methodIndex == 1)
            errordlg('逆衍射逐层测试只支持方法1：反向衍射', '提示');
            set(handles.hStatusText, 'String', '状态：请先切换到方法1');
            drawnow;
            return;
        end

        % 确保预处理全息图可用
        [handles, holoPre] = iEnsurePlaybackHologram(handles, ctx);
        guidata(hFig, handles);

        % 初始化角谱播放状态（准备频域数据）
        try
            state = iInitAngularSpectrumPlayback(holoPre, ctx.zVec, ctx.lambda, ctx.pixel, ctx.useGPU);
        catch ME
            errordlg(sprintf('逆衍射逐层测试初始化失败：\n%s', ME.message), '错误');
            set(handles.hStatusText, 'String', '状态：逆衍射逐层测试初始化失败');
            drawnow;
            return;
        end

        % 创建播放窗口
        hPlayer = figure( ...
            'Name', '逆衍射逐层测试', ...
            'NumberTitle', 'off', ...
            'Color', [0.08 0.09 0.12], ...
            'Position', [180 120 900 720]);

        ax = axes('Parent', hPlayer, ...
            'Units', 'normalized', ...
            'Position', [0.08 0.16 0.84 0.78], ...
            'Color', [0.10 0.11 0.15], ...
            'XColor', [0.9 0.9 0.9], ...
            'YColor', [0.9 0.9 0.9], ...
            'FontName', 'Times New Roman');
        hImg = imagesc(ax, zeros(state.imageSize, 'single'));
        axis(ax, 'image');
        colormap(ax, 'gray');
        colorbar(ax);

        % 信息显示文本
        hInfo = uicontrol('Parent', hPlayer, 'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [0.08 0.05 0.84 0.06], ...
            'BackgroundColor', [0.08 0.09 0.12], ...
            'ForegroundColor', [0.92 0.94 0.98], ...
            'FontName', 'Times New Roman', ...
            'FontSize', 11, ...
            'HorizontalAlignment', 'left', ...
            'String', sprintf('模式：%s    总帧数：%d    关闭此窗口可停止播放。', ...
                              state.mode, state.frameCount));

        set(handles.hStatusText, 'String', '状态：正在播放逆衍射逐层结果...');
        drawnow;

        % 循环播放每一层
        pauseSec = 0.05; % 帧间隔
        stoppedEarly = false;
        for frameIdx = 1:state.frameCount
            % 检查窗口是否被关闭
            if ~isgraphics(hPlayer) || ~isgraphics(hImg)
                stoppedEarly = true;
                break;
            end

            % 计算当前层的衍射场
            frameRaw = iComputeAngularSpectrumPlaybackFrame(state, frameIdx);
            frameDisp = AS_normalizeImage(single(frameRaw)); % 外部函数：归一化显示
            set(hImg, 'CData', frameDisp);
            title(ax, sprintf('逆衍射逐层播放 %d/%d', frameIdx, state.frameCount));
            set(hInfo, 'String', sprintf([ ...
                '模式：%s    当前帧：%d/%d    z = %.2f um    强度范围 = [%.3g, %.3g]'], ...
                state.mode, frameIdx, state.frameCount, state.zVecUm(frameIdx), ...
                min(frameRaw(:)), max(frameRaw(:))));
            drawnow;
            pause(pauseSec);
        end

        handles = guidata(hFig);
        if stoppedEarly
            set(handles.hStatusText, 'String', '状态：逆衍射逐层测试已停止');
        else
            set(handles.hStatusText, 'String', '状态：逆衍射逐层测试已完成');
        end
        drawnow;
    end

    % --- 收集逆衍射播放输入的辅助函数 ---
    function [ctx, errMsg] = iCollectInversePlaybackInputs(handlesIn)
        errMsg = '';
        ctx = struct();

        if isempty(handlesIn.DataImage)
            errMsg = '请先加载全息图。';
            return;
        end

        dataImg = handlesIn.DataImage;
        if ~isempty(handlesIn.BgImage)
            bgImg = handlesIn.BgImage;
            if ~isequal(size(dataImg), size(bgImg))
                bgImg = imresize(bgImg, size(dataImg));
            end
        else
            bgImg = [];
        end

        % 读取并转换参数
        zMin_mm = str2double(get(handlesIn.hZmin, 'String'));
        zMax_mm = str2double(get(handlesIn.hZmax, 'String'));
        dz_um = str2double(get(handlesIn.hDz, 'String'));
        lam_nm = str2double(get(handlesIn.hLam, 'String'));
        pix_um = str2double(get(handlesIn.hPixel, 'String'));

        if any(isnan([zMin_mm, zMax_mm, dz_um, lam_nm, pix_um]))
            errMsg = '请先检查重建参数是否正确。';
            return;
        end

        zMin = zMin_mm * 1e-3;
        zMax = zMax_mm * 1e-3;
        dz = dz_um * 1e-6;
        lambda = lam_nm * 1e-9;
        pixel = pix_um * 1e-6;

        if dz <= 0 || zMax <= zMin
            errMsg = 'z 最小值、最大值和 z 步长无效。';
            return;
        end

        zVec = zMin:dz:zMax;
        if numel(zVec) < 2
            errMsg = '当前 z 范围生成的层数不足 2 层。';
            return;
        end

        % 获取方法名称
        items = get(handlesIn.hMethodPopup, 'String');
        idx = get(handlesIn.hMethodPopup, 'Value');
        if iscell(items)
            methodName = items{idx};
        else
            methodName = items(idx, :);
        end

        % GPU设置
        if isfield(handlesIn, 'hUseGPU') && isgraphics(handlesIn.hUseGPU)
            useGPU = logical(get(handlesIn.hUseGPU, 'Value'));
        else
            useGPU = true;
        end

        ctx.dataImg = dataImg;
        ctx.bgImg = bgImg;
        ctx.zVec = zVec;
        ctx.lambda = lambda;
        ctx.pixel = pixel;
        ctx.methodName = methodName;
        ctx.methodIndex = idx;
        ctx.useGPU = useGPU;
    end

    % --- 确保预处理全息图可用的辅助函数 ---
    function [handlesOut, holoPre] = iEnsurePlaybackHologram(handlesIn, ctx)
        handlesOut = handlesIn;
        hasCachedPreprocess = isfield(handlesIn, 'PreprocessedHologram') && ~isempty(handlesIn.PreprocessedHologram);
        if hasCachedPreprocess
            holoPre = handlesIn.PreprocessedHologram;
            return;
        end

        % 如果没有缓存，则进行预处理
        tPre = tic;
        holoPre = AS_preprocessHologram(ctx.dataImg, ctx.bgImg);
        handlesOut.PreprocessedHologram = holoPre;
        fprintf('[逆衍射逐层测试] 预处理完成，用时 %.2f 秒。\n', toc(tPre));
    end

    % --- 初始化角谱播放状态（频域准备） ---
    function state = iInitAngularSpectrumPlayback(holoPre, zVec, lambda, pixel, requestGPU)
        holo = single(holoPre);
        zVec = single(zVec(:).');
        lambda = single(lambda);
        pixel = single(pixel);

        [ny, nx] = size(holo);
        k = single(2*pi) / lambda; % 波数
        iUnit = complex(single(0), single(1)); % 虚数单位

        useGPU = logical(requestGPU) && iHasUsableGpuLocal();
        state = struct();
        state.imageSize = [ny, nx];
        state.frameCount = numel(zVec);
        state.zVecUm = double(zVec(:)) * 1e6; % 保存微米坐标
        state.mode = 'CPU';

        if useGPU
            try
                g = gpuDevice;
                holoG = gpuArray(holo); % 数据移至GPU
                [FX, FY] = iFrequencyGridLocal(nx, ny, pixel, 'single'); % 生成频域网格
                FX = gpuArray(FX);
                FY = gpuArray(FY);
                lambda2 = lambda * lambda;
                radicand = single(1) - lambda2 .* (FX.^2 + FY.^2);
                state.passband = single(radicand >= 0); % 限制在传播带内（避免倏逝波）
                kz = sqrt(max(radicand, single(0))); % z方向波数
                state.phaseBase = (iUnit * k) .* kz; % 相位基础项
                state.U0f = fft2(holoG); % 对全息图进行傅里叶变换
                state.zVec = gpuArray(zVec);
                wait(g);
                state.mode = 'GPU';
                state.useGPU = true;
                return;
            catch ME
                warning('VolumeGUI_AngularSpectrum:PlaybackGPUFallback', ...
                    '逆衍射逐层测试的 GPU 初始化失败，已回退到 CPU。原因：%s', ...
                    ME.message);
            end
        end

        % CPU 版本
        [FX, FY] = iFrequencyGridLocal(nx, ny, pixel, 'single');
        lambda2 = lambda * lambda;
        radicand = single(1) - lambda2 .* (FX.^2 + FY.^2);
        state.passband = single(radicand >= 0);
        kz = sqrt(max(radicand, single(0)));
        state.phaseBase = (iUnit * k) .* kz;
        state.U0f = fft2(holo);
        state.zVec = zVec;
        state.useGPU = false;
    end

    % --- 计算单帧角谱传播 ---
    function frame = iComputeAngularSpectrumPlaybackFrame(state, frameIdx)
        zNow = state.zVec(frameIdx); % 当前距离
        H = exp(state.phaseBase .* zNow); % 传递函数
        H = H .* state.passband; % 施加通带限制
        Uz = ifft2(state.U0f .* H); % 频域相乘 -> 逆傅里叶变换
        if state.useGPU
            frame = gather(abs(Uz).^2); % 取回CPU并计算强度
        else
            frame = abs(Uz).^2;
        end
    end

    % --- 生成频域坐标网格 ---
    function [FX, FY] = iFrequencyGridLocal(nx, ny, pixel, dtype)
        fx = ifftshift((-floor(nx/2):ceil(nx/2)-1) ./ (double(nx) * double(pixel)));
        fy = ifftshift((-floor(ny/2):ceil(ny/2)-1) ./ (double(ny) * double(pixel)));
        [FX, FY] = meshgrid(cast(fx, dtype), cast(fy, dtype));
    end

    % --- 检查是否有可用GPU ---
    function tf = iHasUsableGpuLocal()
        tf = false;
        try
            if exist('gpuDeviceCount', 'file') == 2
                tf = gpuDeviceCount > 0;
            end
        catch
            tf = false;
        end
    end

    % --- 异步重绘（未在当前UI完全启用） ---
    function onRedrawAsync(~, ~)
        handles = guidata(hFig);
        set(handles.hStatusText, 'String', '状态：已加入3D渲染队列...');
        drawnow;

        if isempty(handles.VolumeData)
            if isfield(handles, 'FastStatsSummary') && ~isempty(handles.FastStatsSummary)
                errordlg('快速统计模式不包含 3D 数据', '提示');
                set(handles.hStatusText, 'String', '状态：快速统计模式不支持3D渲染');
                drawnow;
                return;
            end
            errordlg('请先生成 3D 数据', '提示');
            set(handles.hStatusText, 'String', '状态：重绘失败（无体数据）');
            drawnow;
            return;
        end

        requestAsyncRender('3D 重建结果');
        set(handles.hStatusText, 'String', '状态：已加入3D渲染队列（统计仍可用）');
        drawnow;
    end

    % --- 同步重绘3D ---
    function onRedraw(~, ~)
        handles = guidata(hFig);

        set(handles.hStatusText, 'String', '状态：正在按当前透明度重绘...');
        drawnow;

        if isempty(handles.VolumeData)
            if isfield(handles, 'FastStatsSummary') && ~isempty(handles.FastStatsSummary)
                errordlg('快速统计模式不包含 3D 数据', '提示');
                set(handles.hStatusText, 'String', '状态：快速统计模式不支持3D渲染');
                drawnow;
                return;
            end
            errordlg('请先生成 3D 数据', '提示');
            set(handles.hStatusText, 'String', '状态：重绘失败（无体数据）');
            drawnow;
            return;
        end

        axes(handles.hVolAxes); % 注意：hVolAxes在当前代码片段中未定义，可能是遗留代码
        cla(handles.hVolAxes);
        tDisplay = tic;
        renderCfg = iGetFastRenderConfig(handles);
        renderMeta = show3d(handles.VolumeData, handles.VolumeAlpha, handles.hVolAxes, renderCfg); % 外部函数显示3D
        fprintf('[流程耗时] 3D重绘完成，用时 %.2f 秒。\n', toc(tDisplay));
        fprintf('[3D渲染] 显示尺寸 %dx%dx%d（原始 %dx%dx%d），纹理=%s，步长[%d %d %d]\n', ...
            renderMeta.displaySize(1), renderMeta.displaySize(2), renderMeta.displaySize(3), ...
            renderMeta.sourceSize(1), renderMeta.sourceSize(2), renderMeta.sourceSize(3), ...
            renderMeta.textureMode, renderMeta.stride(1), renderMeta.stride(2), renderMeta.stride(3));
        title(handles.hVolAxes, '3D 重建结果');

        set(handles.hStatusText, 'String', '状态：重绘完成');
        drawnow;
    end

    % --- 请求异步渲染（防抖） ---
    function requestAsyncRender(titleText)
        handles = guidata(hFig);
        if isempty(handles) || isempty(handles.VolumeData)
            return;
        end

        % 如果已有定时器，先停止并删除
        if isfield(handles, 'RenderTimer') && ~isempty(handles.RenderTimer) && isvalid(handles.RenderTimer)
            try
                stop(handles.RenderTimer);
                delete(handles.RenderTimer);
            catch
            end
            handles.RenderTimer = [];
        end

        if ~isfield(handles, 'RenderRequestId') || isempty(handles.RenderRequestId)
            handles.RenderRequestId = 0;
        end
        handles.RenderRequestId = handles.RenderRequestId + 1;
        requestId = handles.RenderRequestId;
        handles.RenderTitle = titleText;
        guidata(hFig, handles);

        % 创建新定时器，延迟0.5秒执行
        t = timer('ExecutionMode', 'singleShot', ...
                  'StartDelay', 0.50, ...
                  'TimerFcn', @(~,~) onAsyncRender(requestId));
        handles = guidata(hFig);
        handles.RenderTimer = t;
        guidata(hFig, handles);
        start(t);
    end

    % --- 异步渲染执行函数 ---
    function onAsyncRender(requestId)
        if ~ishandle(hFig)
            return;
        end
        handles = guidata(hFig);
        % 检查是否是最新的请求
        if isempty(handles) || ~isfield(handles, 'RenderRequestId') || requestId ~= handles.RenderRequestId
            return;
        end
        if isempty(handles.VolumeData)
            return;
        end

        set(handles.hStatusText, 'String', '状态：正在渲染3D（统计仍可用）...');
        drawnow;
        tDisplay = tic;

        try
            axes(handles.hVolAxes);
            cla(handles.hVolAxes);
            renderCfg = iGetFastRenderConfig(handles);
            renderMeta = show3d(handles.VolumeData, handles.VolumeAlpha, handles.hVolAxes, renderCfg);
            if isfield(handles, 'RenderTitle') && ~isempty(handles.RenderTitle)
                title(handles.hVolAxes, handles.RenderTitle);
            end
            if ~handles.ZoomApplied
                try
                    camzoom(handles.hVolAxes, 0.75); % 稍微缩小一点
                catch
                end
                handles.ZoomApplied = true;
            end

            fprintf('[流程耗时] 3D显示完成，用时 %.2f 秒。\n', toc(tDisplay));
            fprintf('[3D渲染] 显示尺寸 %dx%dx%d（原始 %dx%dx%d），纹理=%s，步长[%d %d %d]\n', ...
                renderMeta.displaySize(1), renderMeta.displaySize(2), renderMeta.displaySize(3), ...
                renderMeta.sourceSize(1), renderMeta.sourceSize(2), renderMeta.sourceSize(3), ...
                renderMeta.textureMode, renderMeta.stride(1), renderMeta.stride(2), renderMeta.stride(3));
            set(handles.hStatusText, 'String', '状态：3D渲染完成');
        catch ME
            warning('VolumeGUI_AngularSpectrum:AsyncRenderFailed', ...
                '3D rendering failed: %s', ME.message);
            set(handles.hStatusText, 'String', '状态：3D渲染失败');
        end

        % 清理定时器
        if isfield(handles, 'RenderTimer') && ~isempty(handles.RenderTimer) && isvalid(handles.RenderTimer)
            try
                stop(handles.RenderTimer);
                delete(handles.RenderTimer);
            catch
            end
            handles.RenderTimer = [];
        end

        guidata(hFig, handles);
        drawnow;
    end

    %==================== 3D 显示辅助函数 ====================%
    function renderMeta = show3d(obj_vol, alpham, ax, renderCfg)
        if nargin < 3 || isempty(ax)
            ax = gca;
        end
        if nargin < 4 || isempty(renderCfg)
            renderCfg = iGetFastRenderConfig(guidata(hFig));
        end

        axes(ax);
        cla(ax);

        % 准备显示数据（降采样等）
        [obj_vol, renderMeta] = iPrepareDisplayVolume(obj_vol, renderCfg);
        [Ny, Nx, Nz] = size(obj_vol);
        N = max([Nx, Ny, Nz]);

        set(ax, 'FontSize', 12, ...
            'FontName', 'Times New Roman', ...
            'XTick', 0:floor(N/4):N, ...
            'YTick', 0:floor(N/4):N, ...
            'ZTick', 0:floor(N/4):N, ...
            'XColor', [0.85 0.90 1.00], ...
            'YColor', [0.85 0.90 1.00], ...
            'ZColor', [0.85 0.90 1.00], ...
            'Color', [0.10 0.11 0.15]);

        % 调用外部函数 vol3d 进行体渲染
        vol3d('CData', obj_vol, 'texture', renderMeta.textureMode, 'parent', ax);

        view(ax, -202, 18); % 设置视角
        axis(ax, 'vis3d');
        axis(ax, [0 Nx 0 Ny 0 Nz]);

        colormap(ax, 'hot');
        cm = colormap(ax);
        colormap(ax, 1 - cm); % 反转colormap

        cbar = colorbar('peer', ax, 'Location', 'eastoutside');
        cbar.Color = [1 1 1];
        if isprop(cbar, 'Label')
            cbar.Label.Color = [1 1 1];
        end

        xlabel(ax, 'y', 'FontSize', 16, 'Color', [0.85 0.90 1.00], 'FontAngle', 'italic');
        ylabel(ax, 'x', 'FontSize', 16, 'Color', [0.85 0.90 1.00], 'FontAngle', 'italic');
        zlabel(ax, 'z', 'FontSize', 16, 'Color', [0.85 0.90 1.00], 'FontAngle', 'italic');
        box(ax, 'on');

        % 调整透明度
        fig = ancestor(ax, 'figure');
        alphamap(fig, 'default');
        if alpham ~= 0
            alphamap(fig, 'decrease', alpham);
        end
    end

    %==================== 重建方法与渲染配置 ====================%

    % 获取快速渲染配置
    function renderCfg = iGetFastRenderConfig(handlesIn)
        renderCfg = struct('enabled', true, 'maxDim', 256, 'textureMode', '2D');
        if nargin < 1 || isempty(handlesIn)
            return;
        end
        if isfield(handlesIn, 'FastRenderEnabled') && ~isempty(handlesIn.FastRenderEnabled)
            renderCfg.enabled = logical(handlesIn.FastRenderEnabled);
        end
        if isfield(handlesIn, 'FastRenderMaxDim') && ~isempty(handlesIn.FastRenderMaxDim)
            renderCfg.maxDim = max(64, round(double(handlesIn.FastRenderMaxDim)));
        end
        if isfield(handlesIn, 'FastRenderTexture') && ~isempty(handlesIn.FastRenderTexture)
            textureMode = upper(char(handlesIn.FastRenderTexture));
            if strcmp(textureMode, '2D') || strcmp(textureMode, '3D')
                renderCfg.textureMode = textureMode;
            end
        end
    end

    % 准备用于显示的体数据（降采样）
    function [volDisp, meta] = iPrepareDisplayVolume(volIn, renderCfg)
        srcSize = size(volIn);
        if numel(srcSize) < 3
            srcSize(3) = 1;
        end

        stride = [1, 1, 1];
        textureMode = '3D';
        if renderCfg.enabled
            maxDim = max(64, round(double(renderCfg.maxDim)));
            % 计算降采样步长
            stride(1) = max(1, ceil(srcSize(1) / maxDim));
            stride(2) = max(1, ceil(srcSize(2) / maxDim));
            stride(3) = max(1, ceil(srcSize(3) / maxDim));
            textureMode = renderCfg.textureMode;
        end

        % 执行降采样
        if any(stride > 1)
            volDisp = volIn(1:stride(1):end, 1:stride(2):end, 1:stride(3):end);
        else
            volDisp = volIn;
        end

        if ~isa(volDisp, 'single')
            volDisp = single(volDisp);
        end

        dispSize = size(volDisp);
        if numel(dispSize) < 3
            dispSize(3) = 1;
        end
        meta = struct( ...
            'sourceSize', srcSize(1:3), ...
            'displaySize', dispSize(1:3), ...
            'stride', stride, ...
            'textureMode', textureMode);
    end

    % --- 反向衍射重建（示例函数，代码中未直接调用主流程） ---
    function vol = reconstructInverseDiffraction(dataImg, bgImg, zVec, lambda, pixelSize)
        dataImg = double(dataImg);
        
        if ~isempty(bgImg)
            bgImg = double(bgImg);
            Uo = dataImg - bgImg;
        else
            Uo = dataImg;
        end
       
        lamda = lambda;
        k     = 2*pi/lamda;

        [r, c] = size(Uo);

        Lox = c * pixelSize;
        Loy = r * pixelSize;

        xo = linspace(-Lox/2, Lox/2, c);
        yo = linspace(-Loy/2, Loy/2, r);
        [xo, yo] = meshgrid(xo, yo);

        fa1 = fft2(Uo);

        Nz  = numel(zVec);
        vol = zeros(r, c, Nz);

        for idx = 1:Nz
            zo = -1*zVec(idx);
            % 菲涅尔近似传播
            F0 = exp(1i*k*zo) / (1i*lamda*zo);
            F1 = exp(1i * k/(2*zo) .* (xo.^2 + yo.^2));
            Ff1 = fft2(F1);
            Fuf1 = fa1 .* Ff1;
            Ui = F0 .* fftshift(ifft2(Fuf1));
            Ii = abs(Ui).^2;
            vol(:, :, idx) = Ii;
        end

        maxVal = max(vol(:));
        if maxVal > 0
            vol = vol / maxVal;
        end
    end

    % --- 反卷积重建（占位函数） ---
    function vol = reconstructDeconv(dataImg, bgImg, zVec, lambda, pixelSize) %#ok<INUSD>
        dataImg = double(dataImg);
        if ~isempty(bgImg)
            bgImg = double(bgImg);
            imgNet = dataImg - bgImg;
        else
            imgNet = dataImg;
        end
        imgNet = mat2gray(imgNet);

        h = fspecial('unsharp');
        sharp = imfilter(imgNet, h, 'replicate');

        [ny, ~] = size(imgNet);
        nz = numel(zVec);
        vol = repmat(sharp, [1, 1, nz]); % 简单复制，非真实反卷积
    end

% =========================================================================
    % [新增] 简易粒径统计回调 - 复用主界面数据
    % =========================================================================
    function onSimpleStats(~, ~)
        handles = guidata(hFig); % 获取主界面所有数据
       
        % 1. 数据验证：直接使用主界面共享的 VolumeData
        if isempty(handles.VolumeData)
            errordlg('请先在主界面生成 3D 数据', '无数据');
            return;
        end
        
        % 2. 参数获取：直接使用主界面共享的控件参数
        pix_str = get(handles.hPixel, 'String');
        pix_um = str2double(pix_str);
        if isnan(pix_um), pix_um = 1; end % 默认值防错
       
        % 3. 图像处理 (后台自动完成)
        set(handles.hStatusText, 'String', '状态：正在计算粒子统计...'); drawnow;
        
        % MIP 投影：将 3D 数据压扁为 2D
        volData = handles.VolumeData;
        img2D = max(volData, [], 3);
        
        xyThreshold = str2double(get(handles.hXYThreshold, 'String'));
        if isnan(xyThreshold), xyThreshold = 0.85; end
        levelAbs = min(max(xyThreshold, 0), 1) * max(img2D(:));
        bw = img2D > levelAbs;
        bw = bwareaopen(bw, 5);   % 滤除 < 5像素的噪点
        bw = imfill(bw, 'holes'); % 填充内部空洞
        
        % 4. 统计计算
        stats = regionprops(bw, 'Centroid', 'EquivDiameter');
        
        if isempty(stats)
            msgbox('当前阈值下未检测到粒子', '统计结果');
            set(handles.hStatusText, 'String', '状态：未检测到粒子');
            return;
        end
        
        diams_px = [stats.EquivDiameter];
        diams_um = diams_px * pix_um; % 转换为微米
       
        % 5. 结果展示 (可弹出一个报告窗口)
        hReport = figure('Name', '粒径统计结果报告', ...
            'NumberTitle', 'off', ...
            'Position', [300 300 800 400], ... % 宽一点，并排显示
            'Color', 'white');
            
        % 左图：显示哪里识别成了颗粒(画红圈)
        ax1 = subplot(1, 2, 1);
        imshow(img2D, []); hold on;
        title(ax1, sprintf('识别结果预览'));
        centers = cat(1, stats.Centroid);
        radii = diams_px' / 2;
        viscircles(ax1, centers, radii, 'Color', 'r', 'LineWidth', 0.5);
        
        % 右图：直方图
        ax2 = subplot(1, 2, 2);
        histogram(ax2, diams_um, 20, 'FaceColor', [0.2 0.6 0.8]);
        grid(ax2, 'on');
        title(ax2, '粒径分布');
        xlabel(ax2, ['直径 (' char(181) 'm)']); % 显示 µm
        ylabel(ax2, '粒子数量');
        
        % 底部：文字汇总 (注：原代码此处有语法错误，已修正逻辑)
        infoStr = sprintf(['统计汇总 ', ...
            '粒子总数: %d 个  |   ', ...
            '平均直径: %.2f %sm   |   ', ...
            '最小: %.2f   |   ', ...
            '最大: %.2f'], ...
            length(diams_um), mean(diams_um), char(181), min(diams_um), max(diams_um));
            
        uicontrol('Parent', hReport, 'Style', 'text', ...
            'String', infoStr, ...
            'Units', 'normalized', ...
            'Position', [0 0 1 0.08], ...
            'BackgroundColor', [0.95 0.95 0.95], ...
            'FontSize', 10, 'FontWeight', 'bold');
            
        set(handles.hStatusText, 'String', '状态：统计报告已生成');
    end

    % --- V2版本统计主函数（当前使用） ---
    function onSimpleStatsV2(~, ~)
        handles = guidata(hFig);
        [handles, statsResult, errMsg] = iBuildStatsResult(handles);
        if ~isempty(errMsg)
            errordlg(errMsg, '统计结果');
            if contains(errMsg, '未检测到粒子')
                set(handles.hStatusText, 'String', '状态：未检测到粒子');
            else
                set(handles.hStatusText, 'String', '状态：统计结果生成失败');
            end
            handles = iRefreshUiState(handles);
            guidata(hFig, handles);
            return;
        end

        handles.StatsResult = statsResult;
        handles.ParticleCoords3D = statsResult.coordsPx;
        handles.ParticleAxialCurves = {};
        handles = iRefreshUiState(handles);
        iUpdateStatsDisplay(handles);
        guidata(hFig, handles);

        iShowParticleOverview(statsResult); % 弹出概览窗口
        set(handles.hStatusText, 'String', sprintf('状态：已完成统计，粒子数%d', size(statsResult.coordsPx, 1)));
    end

    % --- 导出结果回调 ---
    function onExportResults(~, ~)
        handles = guidata(hFig);
        if ~isfield(handles, 'StatsResult') || isempty(handles.StatsResult)
            errordlg('当前没有可导出的统计结果', '导出结果');
            return;
        end

        [file, path, filterIdx] = uiputfile( ...
            {'*.xlsx', 'Excel 文件 (*.xlsx)'; '*.csv', 'CSV 文件 (*.csv)'}, ...
            '导出粒子统计结果');
        if isequal(file, 0)
            return;
        end

        targetPath = fullfile(path, file);
        try
            iExportStatsResult(handles.StatsResult, targetPath, filterIdx);
            set(handles.hStatusText, 'String', '状态：导出完成');
        catch ME
            errordlg(sprintf('导出失败：\n%s', ME.message), '导出结果');
            set(handles.hStatusText, 'String', '状态：导出失败');
        end
        drawnow;
    end

    % --- 方法切换回调 ---
    function onMethodChanged(~, ~)
        handles = guidata(hFig);
        handles = iRefreshUiState(handles); % 刷新控件状态
        guidata(hFig, handles);
    end

    % --- 刷新UI控件状态（启用/禁用） ---
    function handlesOut = iRefreshUiState(handlesIn)
        handlesOut = handlesIn;

        % 方法1保留逆衍射测试入口，导出按钮由统计结果控制。
        if isfield(handlesIn, 'hExportBtn') && isgraphics(handlesIn.hExportBtn)
            exportState = 'off';
            if isfield(handlesIn, 'StatsResult') && ~isempty(handlesIn.StatsResult)
                exportState = 'on';
            end
            set(handlesIn.hExportBtn, 'Enable', exportState);
        end

        if isgraphics(hInverseTest)
            set(hInverseTest, 'Enable', 'on');
        end
    end
    % --- 更新统计显示（直方图、表格、摘要） ---
    function iUpdateStatsDisplay(handlesIn)
        axes(handlesIn.hHistAxes);
        cla(handlesIn.hHistAxes);

        % 如果没有结果，清空显示
        if ~isfield(handlesIn, 'StatsResult') || isempty(handlesIn.StatsResult)
            title(handlesIn.hHistAxes, '粒径分布', 'Color', textColor);
            xlabel(handlesIn.hHistAxes, ['直径 (' muChar 'm)'], 'Color', textColor);
            ylabel(handlesIn.hHistAxes, '粒子数量', 'Color', textColor);
            grid(handlesIn.hHistAxes, 'on');
            set(handlesIn.hResultTable, 'Data', cell(0, 5));
            set(handlesIn.hSummaryText, 'String', '尚未生成统计结果');
            return;
        end

        statsResult = handlesIn.StatsResult;
        % 画直方图
        histogram(handlesIn.hHistAxes, statsResult.diamsUm, 20, 'FaceColor', [0.2 0.6 0.8]);
        set(handlesIn.hHistAxes, 'Color', [0.10 0.11 0.15], ...
            'XColor', textColor, 'YColor', textColor, 'FontName', 'Times New Roman');
        grid(handlesIn.hHistAxes, 'on');
        title(handlesIn.hHistAxes, '粒径分布', 'Color', textColor);
        xlabel(handlesIn.hHistAxes, ['直径 (' muChar 'm)'], 'Color', textColor);
        ylabel(handlesIn.hHistAxes, '粒子数量', 'Color', textColor);

        % 更新摘要文本
        summaryLines = { ...
            sprintf('方法: %s', statsResult.methodName), ...
            sprintf('粒子总数: %d', statsResult.count), ...
            sprintf('平均粒径: %.2f %cm', statsResult.meanDiamUm, muChar), ...
            sprintf('MVD: %.2f %cm', statsResult.mvdUm, muChar), ...
            sprintf('D[3,2]: %.2f %cm', statsResult.meanEffDiamUm, muChar), ...
            sprintf('浓度: %.3g 个/mL', statsResult.numConcPerMl), ...
            sprintf('z范围: %.0f ~ %.0f %cm', statsResult.zRangeUm(1), statsResult.zRangeUm(2), muChar), ...
            sprintf('dz: %.0f %cm', statsResult.zStepUm, muChar)};
        set(handlesIn.hSummaryText, 'String', summaryLines);

        % 更新表格
        tableData = num2cell([(1:size(statsResult.coordsUm, 1)).', statsResult.coordsUm, statsResult.diamsUm(:)]);
        set(handlesIn.hResultTable, 'Data', tableData);
    end

    % --- 构建统计结果的核心逻辑 ---
    function [handlesOut, statsResult, errMsg] = iBuildStatsResult(handlesIn)
        handlesOut = handlesIn;
        statsResult = [];
        errMsg = '';

        fastStats = [];
        if isfield(handlesIn, 'FastStatsSummary') && ~isempty(handlesIn.FastStatsSummary)
            fastStats = handlesIn.FastStatsSummary;
        end

        if isempty(fastStats)
            errMsg = '请先完成方法1重建。';
            return;
        end

        pix_um = str2double(get(handlesIn.hPixel, 'String'));
        if isnan(pix_um)
            pix_um = 1;
        end

        coordsPx = fastStats.coords3D;
        diamsPx = fastStats.diamsPx(:);
        diamsUm = diamsPx * pix_um;
        zVecUm = fastStats.zVecUm(:);
        reconRows = fastStats.volumeSize(1);
        reconCols = fastStats.volumeSize(2);
        if isempty(coordsPx) || isempty(diamsUm)
            errMsg = '当前阈值下未检测到粒子，请检查重建结果。';
            return;
        end

        zStepUm = 1;
        if numel(zVecUm) > 1
            zStepUm = median(abs(diff(zVecUm)));
        end
        if ~(isfinite(zStepUm) && zStepUm > 0)
            zStepUm = 1;
        end
        zSpanUm = max(zStepUm, abs(zVecUm(end) - zVecUm(1)) + zStepUm);
        reconVolumeUm3 = reconRows * reconCols * (pix_um ^ 2) * zSpanUm;
        numConcPerMl = size(coordsPx, 1) * 1e12 / max(reconVolumeUm3, eps); % 个/mL

        diamVecUm = diamsUm(:);
        volWeights = diamVecUm .^ 3;
        if any(volWeights > 0)
            [diamSorted, sortIdx] = sort(diamVecUm);
            cumFrac = cumsum(volWeights(sortIdx)) / sum(volWeights);
            idx50 = find(cumFrac >= 0.5, 1, 'first'); % 体积中位径
            if isempty(idx50)
                mvdUm = diamSorted(end);
            elseif idx50 == 1
                mvdUm = diamSorted(1);
            else
                prevFrac = cumFrac(idx50 - 1);
                nextFrac = cumFrac(idx50);
                prevDiam = diamSorted(idx50 - 1);
                nextDiam = diamSorted(idx50);
                if nextFrac <= prevFrac
                    interpT = 0;
                else
                    interpT = (0.5 - prevFrac) / (nextFrac - prevFrac);
                end
                mvdUm = prevDiam + interpT * (nextDiam - prevDiam);
            end
        else
            mvdUm = median(diamVecUm);
        end

        denomD32 = sum(diamVecUm .^ 2);
        if denomD32 > 0
            meanEffDiamUm = sum(diamVecUm .^ 3) / denomD32;
        else
            meanEffDiamUm = mean(diamVecUm);
        end

        coordsUm = [coordsPx(:, 1) * pix_um, coordsPx(:, 2) * pix_um, coordsPx(:, 3)];
        methodName = '方法1：反向衍射';
        mipImage = fastStats.img2D;
        roiBoxesAll = fastStats.candidateRoiBoxes;
        validMaskAll = logical(fastStats.validMask(:));
        candidateCoordsPx = fastStats.candidateCoords3D;

        statsResult = struct( ...
            'methodName', methodName, ...
            'coordsPx', coordsPx, ...
            'coordsUm', coordsUm, ...
            'diamsUm', diamsUm(:), ...
            'diamsPx', diamsPx(:), ...
            'pixUm', pix_um, ...
            'zVecUm', zVecUm(:), ...
            'zStepUm', zStepUm, ...
            'zRangeUm', [zVecUm(1), zVecUm(end)], ...
            'count', size(coordsPx, 1), ...
            'meanDiamUm', mean(diamVecUm), ...
            'mvdUm', mvdUm, ...
            'meanEffDiamUm', meanEffDiamUm, ...
            'numConcPerMl', numConcPerMl, ...
            'mipImage', mipImage, ...
            'candidateCoordsPx', candidateCoordsPx, ...
            'candidateRoiBoxes', roiBoxesAll, ...
            'candidateValidMask', validMaskAll, ...
            'reconRows', reconRows, ...
            'reconCols', reconCols, ...
            'paramSnapshot', struct('pix_um', pix_um, 'z_min_um', zVecUm(1), 'z_max_um', zVecUm(end), 'dz_um', zStepUm));
    end

    % --- 显示粒子概览窗口（MIP+3D散点） ---
    function hOverview = iShowParticleOverview(statsResult)
        if isempty(statsResult) || ~isfield(statsResult, 'mipImage') || isempty(statsResult.mipImage)
            hOverview = [];
            return;
        end

        hOverview = figure( ...
            'Name', '粒子识别总览', ...
            'NumberTitle', 'off', ...
            'Color', [0.05 0.07 0.12], ...
            'Position', [220 120 1360 620]);

        % 左轴：MIP 图与粒子标注
        axMip = axes('Parent', hOverview, ...
            'Units', 'normalized', ...
            'Position', [0.05 0.12 0.42 0.80], ...
            'Color', [0.10 0.11 0.15], ...
            'XColor', [0.85 0.90 1.00], ...
            'YColor', [0.85 0.90 1.00], ...
            'FontName', 'Times New Roman');
        imagesc(axMip, statsResult.mipImage);
        axis(axMip, 'image');
        colormap(axMip, 'gray');
        title(axMip, 'MIP 候选粒子与有效粒子识别', 'Color', [0.90 0.95 1.00]);
        hold(axMip, 'on');

        roiBoxes = statsResult.candidateRoiBoxes;
        validMask = logical(statsResult.candidateValidMask(:));

        % --- 画候选 ROI 框 (candidateRoiBoxes)，而不是等效直径方框 ---
        % 仿真模式下 ROI = BoundingBox + margin，实验模式下 ROI = 自适应扩展。
        % 用 candidateRoiBoxes 才能看到真实候选区域是否完整框住粒子。
        validCount = 0;
        for n = 1:size(roiBoxes, 1)
            if validMask(n)
                boxColor = [0.20 0.95 0.55];  % green: valid
                validCount = validCount + 1;
            else
                boxColor = [1.00 0.65 0.20];  % orange: rejected
            end
            rectangle(axMip, 'Position', roiBoxes(n, :), ...
                'EdgeColor', boxColor, 'LineWidth', 1.5);

            % Label at top-left corner of each box
            text(axMip, roiBoxes(n, 1) + 2, roiBoxes(n, 2) - 4, sprintf('%d', n), ...
                'Color', boxColor, 'FontSize', 10, 'FontWeight', 'bold', 'Parent', axMip);
        end

        validCount = sum(validMask);
        candidateCount = size(roiBoxes, 1);
        text(axMip, 0.02, 0.98, sprintf('候选粒子:%d    有效粒子: %d', candidateCount, validCount), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'Color', [0.95 0.97 1.00], 'FontSize', 11, 'FontWeight', 'bold', 'Parent', axMip);

        % 右轴：3D 散点/球体
        ax3D = axes('Parent', hOverview, ...
            'Units', 'normalized', ...
            'Position', [0.54 0.12 0.41 0.80], ...
            'Color', [0.10 0.11 0.15], ...
            'XColor', [0.85 0.90 1.00], ...
            'YColor', [0.85 0.90 1.00], ...
            'ZColor', [0.85 0.90 1.00], ...
            'FontName', 'Times New Roman');
        AS_showParticleSpheres(statsResult.coordsPx, statsResult.diamsUm, statsResult.pixUm, ax3D); % 外部函数
        title(ax3D, '稀疏粒子场', 'Color', [0.90 0.95 1.00]);
    end

    % --- 导出统计结果到文件 ---
    function iExportStatsResult(statsResult, targetPath, filterIdx)
        % 构建粒子表格
        particleTbl = table((1:statsResult.count).', statsResult.coordsUm(:, 1), statsResult.coordsUm(:, 2), ...
            statsResult.coordsUm(:, 3), statsResult.diamsUm(:), ...
            'VariableNames', {'ID', 'x_um', 'y_um', 'z_um', 'd_um'});
        % 构建摘要表格
        summaryRows = { ...
            'method', statsResult.methodName; ...
            'particle_count', statsResult.count; ...
            'mean_diameter_um', statsResult.meanDiamUm; ...
            'mvd_um', statsResult.mvdUm; ...
            'd32_um', statsResult.meanEffDiamUm; ...
            'concentration_per_mL', statsResult.numConcPerMl; ...
            'pix_um', statsResult.paramSnapshot.pix_um; ...
            'z_min_um', statsResult.paramSnapshot.z_min_um; ...
            'z_max_um', statsResult.paramSnapshot.z_max_um; ...
            'dz_um', statsResult.paramSnapshot.dz_um};
        summaryTbl = cell2table(summaryRows, 'VariableNames', {'Metric', 'Value'});

        [folder, baseName, ext] = fileparts(targetPath);
        if isempty(ext)
            if filterIdx == 1
                ext = '.xlsx';
            else
                ext = '.csv';
            end
        end

        % 写入 Excel
        if strcmpi(ext, '.xlsx') || filterIdx == 1
            excelPath = fullfile(folder, [baseName, '.xlsx']);
            writetable(particleTbl, excelPath, 'Sheet', 'Particles');
            writetable(summaryTbl, excelPath, 'Sheet', 'Summary');
        else
            % 写入 CSV
            particlePath = fullfile(folder, [baseName, '_particles.csv']);
            summaryPath = fullfile(folder, [baseName, '_summary.csv']);
            writetable(particleTbl, particlePath);
            writetable(summaryTbl, summaryPath);
        end
    end

    % --- 重绘V2（未完全使用） ---
    function onRedrawV2(~, ~)
        handles = guidata(hFig);

        set(handles.hStatusText, 'String', '状态：正在按当前透明度重绘...');
        drawnow;

        if isempty(handles.VolumeData)
            if isfield(handles, 'FastStatsSummary') && ~isempty(handles.FastStatsSummary)
                errordlg('快速统计模式不包含 3D 数据', '提示');
                set(handles.hStatusText, 'String', '状态：快速统计模式不支持3D渲染');
                drawnow;
                return;
            end
            errordlg('请先生成 3D 数据', '提示');
            set(handles.hStatusText, 'String', '状态：重绘失败（无体数据）');
            drawnow;
            return;
        end

        onRedraw([], []);
        set(handles.hStatusText, 'String', '状态：重绘完成');
        drawnow;
    end
end
