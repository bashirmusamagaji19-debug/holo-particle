function VolumeGUI_Classic
% 3 列经典 GUI：
% 左：参数设置
% 中：数据图 / 背景图
% 右：3D 体绘制（使用 show3d + alphamap）
    addpath('./function/');

    %==================== figure ====================%
    hFig = figure( ...
        'Name', '3D Digital Hologram Reconstruction', ...
        'NumberTitle', 'off', ...
        'Color', [0.05 0.07 0.12], ...
        'Position', [100 100 1400 600]);

    % 顶部白色标题
    uicontrol('Parent', hFig, 'Style', 'text', ...
        'String', '3D Digital Hologram Reconstruction', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.92 0.6 0.06], ...
        'ForegroundColor', [0.90 0.95 1.00], ...
        'BackgroundColor', [0.05 0.07 0.12], ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left');

    % 简单“主题色”
    primaryColor   = [0.00 0.65 0.70];
    secondaryColor = [0.16 0.32 0.63];
    accentColor    = [0.30 0.50 0.80];
    textColor      = [0.85 0.90 1.00];
    panelBg        = [0.09 0.10 0.14];

    % 一些特殊字符
    muChar  = char(181);   % µ
    lamChar = char(955);   % λ

    % 结构保存数据
    handles = struct();
    handles.DataImage     = [];
    handles.BgImage       = [];
    handles.VolumeData    = [];
    handles.VolumeZVecUm  = [];
    handles.ParticleCoords3D = [];
    handles.ParticleAxialCurves = {};
    handles.VolumeAlpha   = 0.1;   % alphamap('decrease', alpham) 的 alpham
    handles.ZoomApplied   = false; % 是否已经调用过 camzoom

  %==================== figure ====================%
    hFig = figure( ...
        'Name', '3D Digital Hologram Reconstruction', ...
        'NumberTitle', 'off', ...
        'Color', [0.05 0.07 0.12], ...
        'Position', [100 100 1400 600]);

    % --- [新增] 顶部菜单栏 ---
    hMenu = uimenu(hFig, 'Label', '分析 (Analysis)');
    uimenu(hMenu, ...
        'Label', '生成粒径统计报告 (Generate Particle Report)', ...
        'Callback', @onSimpleStatsV2);
    % ------------------------
    


    %==================== 左列：参数 panel ====================%
    leftPanel = uipanel('Parent', hFig, ...
        'Title', ' Reconstruction Parameters ', ...
        'Units', 'normalized', ...
        'Position', [0.01 0.10 0.18 0.8], ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', textColor, ...
        'BorderType', 'line', ...
        'HighlightColor', [0.4 0.4 0.4], ...
        'ShadowColor', [0 0 0], ...
        'FontSize', 12, 'FontWeight', 'bold');

    % z_min (mm)
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'z-min (mm):', ...
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
        'String', '0');

    % z_max (mm)
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'z-max (mm):', ...
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
        'String', '6');

    % dz (µm)
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', ['dz (' muChar 'm):'], ...
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
        'String', '50');   % 你原 CS 代码里 deltaZ=50

    % λ (nm)
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', [lamChar ' (nm):'], ...
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
        'String', '632.8');  % 例：660 nm≈0.66 µm

    % PixelSize (µm) -> detector_size
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', ['Pixel Size (' muChar 'm):'], ...
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
        'String', '2');    % 你原 CS 里 detector_size=2

    % 重建方法
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
        'String', { ...
            '方法1：逆向衍射', ...
            '方法2：压缩感知'});

    % ---- CS 参数：pad_size ----
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'Padsize:', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.34 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hPadSize = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.34 0.45 0.06], ...
        'String', '0');  % 你原代码 pad_size=0

    % tau
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'tau:', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.27 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hTau = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.27 0.45 0.06], ...
        'String', '0.01');

    % piter
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'piter:', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.20 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hPiter = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.20 0.45 0.06], ...
        'String', '2');

    % tolA
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'tolA:', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.13 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hTolA = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.13 0.45 0.06], ...
        'String', '1e-5');

    % iterations
    uicontrol('Parent', leftPanel, 'Style', 'text', ...
        'String', 'iterations:', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.06 0.40 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', panelBg, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'left');
    handles.hIter = uicontrol('Parent', leftPanel, 'Style', 'edit', ...
        'Units', 'normalized', ...
        'Position', [0.48 0.06 0.45 0.06], ...
        'String', '100');

    %==================== 中列：数据/背景 ====================%
    set(handles.hDz, 'String', '30');

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

    handles.hDataAxes = axes('Parent', midPanel, ...
        'Units', 'normalized', ...
        'Position', [0.10 0.60 0.85 0.34]);
    title(handles.hDataAxes, '数据图 1');
    set(handles.hDataAxes, 'Color', [0.10 0.11 0.15], ...
        'XColor', textColor, 'YColor', textColor, 'FontName', 'Times New Roman');

    hLoadData = uicontrol('Parent', midPanel, 'Style', 'pushbutton', ...
        'String', '加载数据图1', ...
        'Units', 'normalized', ...
        'Position', [0.25 0.48 0.50 0.06], ...
        'BackgroundColor', secondaryColor, ...
        'ForegroundColor', [0.95 0.97 1.00], ...
        'FontWeight', 'bold');

    handles.hBgAxes = axes('Parent', midPanel, ...
        'Units', 'normalized', ...
        'Position', [0.10 0.12 0.85 0.34]);
    title(handles.hBgAxes, '背景图 2');
    set(handles.hBgAxes, 'Color', [0.10 0.11 0.15], ...
        'XColor', textColor, 'YColor', textColor, 'FontName', 'Times New Roman');

    hLoadBg = uicontrol('Parent', midPanel, 'Style', 'pushbutton', ...
        'String', '加载背景图2', ...
        'Units', 'normalized', ...
        'Position', [0.25 0.01 0.50 0.06], ...
        'BackgroundColor', secondaryColor, ...
        'ForegroundColor', [0.95 0.97 1.00], ...
        'FontWeight', 'bold');

    %==================== 右列：3D 显示 ====================%
    rightPanel = uipanel('Parent', hFig, ...
        'Title', ' 3D Reconstruction Volume ', ...
        'Units', 'normalized', ...
        'Position', [0.50 0.10 0.49 0.8], ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'ForegroundColor', textColor, ...
        'BorderType', 'line', ...
        'HighlightColor', [0.4 0.4 0.4], ...
        'ShadowColor', [0 0 0], ...
        'FontSize', 12, 'FontWeight', 'bold');

    % 3D 轴：右侧占大部分区域，左边留出位置给按钮
    handles.hVolAxes = axes('Parent', rightPanel, ...
        'Units', 'normalized', ...
        'Position', [0.38  0.15  0.60  0.80]);
%     title(handles.hVolAxes, '3D 重建结果');
    xlabel(handles.hVolAxes, 'y', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
    ylabel(handles.hVolAxes, 'x', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
    zlabel(handles.hVolAxes, 'z', 'FontName', 'Times New Roman', 'FontAngle', 'italic');
    set(handles.hVolAxes, 'Color', [0.10 0.11 0.15], ...
        'XColor', textColor, 'YColor', textColor, 'ZColor', textColor, ...
        'FontName', 'Times New Roman');

    % 左侧竖排：Alpha + 两个按钮
    uicontrol('Parent', rightPanel, 'Style', 'text', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'String', 'Alpha (decrease):', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.78 0.22 0.06], ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'HorizontalAlignment', 'left');

    handles.hAlphaSlider = uicontrol('Parent', rightPanel, 'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.73 0.22 0.04], ...
        'Min', 0, 'Max', 1, 'Value', handles.VolumeAlpha, ...
        'BackgroundColor', [0.15 0.17 0.22]);

    handles.hAlphaText = uicontrol('Parent', rightPanel, 'Style', 'text', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'Units', 'normalized', ...
        'Position', [0.03 0.69 0.22 0.03], ...
        'String', sprintf('Alpha = %.2f', handles.VolumeAlpha), ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'HorizontalAlignment', 'left');

    hRedraw = uicontrol('Parent', rightPanel, 'Style', 'pushbutton', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.55 0.22 0.08], ...
        'String', '按当前透明度重绘', ...
        'BackgroundColor', accentColor, ...
        'ForegroundColor', [0.05 0.07 0.10], ...
        'FontWeight', 'bold');

    hReconstruct = uicontrol('Parent', rightPanel, 'Style', 'pushbutton', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.42 0.22 0.08], ...
        'String', '计算并绘制3D图', ...
        'BackgroundColor', primaryColor, ...
        'ForegroundColor', [0.05 0.07 0.10], ...
        'FontWeight', 'bold');

    % 状态显示文本框（新增）
    handles.hStatusText = uicontrol('Parent', rightPanel, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.03 0.30 0.35 0.08], ...
        'String', 'Status: Idle', ...
        'ForegroundColor', textColor, ...
        'BackgroundColor', [0.07 0.08 0.11], ...
        'FontName', 'Times New Roman', ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'left');

    guidata(hFig, handles);

    %==================== 绑定回调 ====================%
    set(hLoadData,    'Callback', @onLoadData);
    set(hLoadBg,      'Callback', @onLoadBg);
    set(hReconstruct, 'Callback', @onReconstruct);
    set(hRedraw,      'Callback', @onRedraw);
    set(handles.hAlphaSlider, 'Callback', @onAlphaChanged);

    %==================== 嵌套函数 ====================%

    function onLoadData(~, ~)
        handles = guidata(hFig);
        [file, path] = uigetfile( ...
            {'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Image Files'}, ...
            '选择数据图片(图片1)');
        if isequal(file, 0), return; end

        img = imread(fullfile(path, file));
        
        if ndims(img) == 3, img = rgb2gray(img); end
        img = double(img);

        handles.DataImage  = img;
        handles.VolumeData = [];
        handles.VolumeZVecUm = [];
        handles.ParticleCoords3D = [];
        handles.ParticleAxialCurves = {};
        handles.ZoomApplied = false;  % 换数据后，下次重建重新 zoom 一次

        axes(handles.hDataAxes);
        cla(handles.hDataAxes);

        imagesc(handles.hDataAxes, img);
        axis(handles.hDataAxes, 'image');
        colormap(handles.hDataAxes, 'gray');

        % colorbar + 字体
        cbar = colorbar(handles.hDataAxes);
        cbar.FontName = 'Times New Roman';
        cbar.FontSize = 9;
        cbar.Color    = [1 1 1];

        set(handles.hDataAxes, ...
            'FontName', 'Times New Roman', ...
            'FontSize', 9, ...
            'XColor', [1 1 1], ...
            'YColor', [1 1 1]);

        % 状态更新
        set(handles.hStatusText, 'String', 'Status: Hologram loaded');
        drawnow;

        guidata(hFig, handles);
    end

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
        handles.ParticleCoords3D = [];
        handles.ParticleAxialCurves = {};

        axes(handles.hBgAxes);
        cla(handles.hBgAxes);
        imagesc(handles.hBgAxes, img);
        axis(handles.hBgAxes, 'image');
        colormap(handles.hBgAxes, 'gray');
        colorbar(handles.hBgAxes);
        title(handles.hBgAxes, '背景图 2');

        % 状态更新
        set(handles.hStatusText, 'String', 'Status: Background Image loaded');
        drawnow;

        guidata(hFig, handles);
    end

    function onAlphaChanged(src, ~)
        handles = guidata(hFig);
        val = get(src, 'Value');  % 0~1
        handles.VolumeAlpha = val;
        set(handles.hAlphaText, 'String', sprintf('Alpha = %.2f', val));
        guidata(hFig, handles);
    end

    function onReconstruct(~, ~)
        handles = guidata(hFig);

        % 状态：开始重建
        set(handles.hStatusText, 'String', 'Status: Reconstructing 3D volume...');
        drawnow;

        if isempty(handles.DataImage)
            errordlg('请先加载数据图1', '提示');
            set(handles.hStatusText, 'String', 'Status: Reconstruction failed (no data image loaded)');
            drawnow;
            return;
        end

        dataImg = handles.DataImage;

        if ~isempty(handles.BgImage)
            bgImg = handles.BgImage;
            if ~isequal(size(dataImg), size(bgImg))
                bgImg = imresize(bgImg, size(dataImg));
            end
        else
            bgImg = [];
        end

        % 读取基础参数
        zMin_mm = str2double(get(handles.hZmin,  'String'));
        zMax_mm = str2double(get(handles.hZmax,  'String'));
        dz_um   = str2double(get(handles.hDz,    'String'));
        lam_nm  = str2double(get(handles.hLam,   'String'));
        pix_um  = str2double(get(handles.hPixel, 'String'));

        % CS 参数
        pad_size  = str2double(get(handles.hPadSize, 'String'));
        tau       = str2double(get(handles.hTau,     'String'));
        piter     = str2double(get(handles.hPiter,   'String'));
        tolA      = str2double(get(handles.hTolA,    'String'));
        iterations= str2double(get(handles.hIter,    'String'));

        if any(isnan([zMin_mm, zMax_mm, dz_um, lam_nm, pix_um, ...
                      pad_size, tau, piter, tolA, iterations]))
            errordlg('请检查所有参数是否为数字', '参数错误');
            set(handles.hStatusText, 'String', 'Status: Reconstruction failed (parameter error)');
            drawnow;
            return;
        end

        zMin   = zMin_mm * 1e-3;   % mm -> m
        zMax   = zMax_mm * 1e-3;   % mm -> m
        dz     = dz_um   * 1e-6;   % µm -> m
        lambda = lam_nm  * 1e-9;   % nm -> m
        pixel  = pix_um  * 1e-6;   % µm -> m

        if dz <= 0 || zMax <= zMin
            errordlg('请检查 z_min, z_max 和 dz：z_max > z_min 且 dz > 0', '参数错误');
            set(handles.hStatusText, 'String', 'Status: Reconstruction failed (invalid z parameters)');
            drawnow;
            return;
        end

        zVec = zMin:dz:zMax;
        if numel(zVec) < 2
            errordlg('z 范围太小或步长太大，导致层数不足', '参数错误');
            set(handles.hStatusText, 'String', 'Status: Reconstruction failed (insufficient z-layers)');
            drawnow;
            return;
        end

        items = get(handles.hMethodPopup, 'String');
        idx   = get(handles.hMethodPopup, 'Value');
        methodName = items{idx};

        switch methodName
            case '方法1：逆向衍射'
                vol = reconstructInverseDiffraction(dataImg, bgImg, zVec, lambda, pixel);

            case '方法2：压缩感知'
                vol = reconstructCS(dataImg, bgImg, ...
                                    zMin_mm, zMax_mm, dz_um, ...
                                    lam_nm, pix_um, ...
                                    pad_size, tau, piter, tolA, iterations);

            otherwise
                errordlg('未知方法选择', '错误');
                set(handles.hStatusText, 'String', 'Status: Reconstruction failed (unknown method)');
                drawnow;
                return;
        end

        if ndims(vol) ~= 3
            errordlg('重建函数没有返回 3D 体数据（尺寸不为 N×M×K）', '错误');
            set(handles.hStatusText, 'String', 'Status: Reconstruction failed (result is not a 3D volume)');
            drawnow;
            return;
        end

        handles.VolumeData = vol;
        handles.VolumeZVecUm = zMin_mm * 1000 + (0:size(vol, 3)-1) * dz_um;
        handles.ParticleCoords3D = [];
        handles.ParticleAxialCurves = {};

        axes(handles.hVolAxes);
        cla(handles.hVolAxes);
        show3d(vol, handles.VolumeAlpha, handles.hVolAxes);
        title(handles.hVolAxes, methodName);

        % 只在第一次重建时 zoom 一次
        if ~handles.ZoomApplied
            try
                camzoom(handles.hVolAxes, 0.75);
            catch
            end
            handles.ZoomApplied = true;
        end

        % 状态：重建完成
        set(handles.hStatusText, 'String', 'Status: 3D reconstruction finished');
        drawnow;

        guidata(hFig, handles);
    end

    function onRedraw(~, ~)
        handles = guidata(hFig);

        % 状态：开始重绘
        set(handles.hStatusText, 'String', 'Status: Redrawing with current alpha...');
        drawnow;

        if isempty(handles.VolumeData)
            errordlg('请先点击“计算并绘制3D图”得到体数据', '提示');
            set(handles.hStatusText, 'String', 'Status: Redraw failed (no volume data)');
            drawnow;
            return;
        end

        axes(handles.hVolAxes);
        cla(handles.hVolAxes);
        show3d(handles.VolumeData, handles.VolumeAlpha, handles.hVolAxes);
        title(handles.hVolAxes, '3D 重建结果');

        % 状态：重绘完成
        set(handles.hStatusText, 'String', 'Status: Redraw finished');
        drawnow;
    end

    %==================== show3d ====================%
    function show3d(obj_vol, alpham, ax)
        if nargin < 3 || isempty(ax)
            ax = gca;
        end

        axes(ax);
        cla(ax);

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

        vol3d('CData', obj_vol, 'texture', '3D', 'parent', ax);

        view(ax, -202, 18);
        axis(ax, 'vis3d');
        axis(ax, [0 Nx 0 Ny 0 Nz]);

        colormap(ax, 'hot');
        cm = colormap(ax);
        colormap(ax, 1 - cm);

        cbar = colorbar('peer', ax, 'Location', 'eastoutside');
        cbar.Color = [1 1 1];
        if isprop(cbar, 'Label')
            cbar.Label.Color = [1 1 1];
        end

        xlabel(ax, 'y', 'FontSize', 16, 'Color', [0.85 0.90 1.00], 'FontAngle', 'italic');
        ylabel(ax, 'x', 'FontSize', 16, 'Color', [0.85 0.90 1.00], 'FontAngle', 'italic');
        zlabel(ax, 'z', 'FontSize', 16, 'Color', [0.85 0.90 1.00], 'FontAngle', 'italic');
        box(ax, 'on');

        % 透明度
        fig = ancestor(ax, 'figure');
        alphamap(fig, 'default');
        if alpham ~= 0
            alphamap(fig, 'decrease', alpham);
        end
    end

    %==================== 三种重建方法 ====================%

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
            zo = -1*zVec(idx)
            lamda
            pixelSize
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
        % %{--- 新增：硬阈值过滤 ---
        %  threshold = 0.7; % 过滤掉低于最大强度 15% 的噪声
        %   vol(vol < threshold) = 0; 
          
    end

    % ==== 压缩感知 + TwIST ====
    function vol = reconstructCS(dataImg, bgImg, ...
                                 zMin_mm, zMax_mm, dz_um, ...
                                 lam_nm, pix_um, ...
                                 pad_size, tau, piter, tolA, iterations)

        % 确保 Functions 在路径上（MyMakingPhase3D / TwIST 等）
        if exist('./Functions', 'dir')
            addpath('./Functions');
        end

        % hologram：dataImg 减背景
        holo = double(dataImg);
        if ~isempty(bgImg)
            bgImg = double(bgImg);
            if ~isequal(size(holo), size(bgImg))
                bgImg = imresize(bgImg, size(holo));
            end
            holo = holo - bgImg;
        end

        % pad_size 取非负整数
        pad_size = max(0, round(pad_size));

        pixel_num    = size(holo, 1); % 假设正方探测器
        detector_size= pix_um;        % µm
        lambda_um    = lam_nm / 1000; % nm -> µm
        deltaZ       = dz_um;         % µm
        offsetZ      = zMin_mm * 1000;% mm -> µm

        totalDepth_um = max(0, (zMax_mm - zMin_mm) * 1000);
        nz = max(1, floor(totalDepth_um / deltaZ) + 1);

        g = holo;

        shrinkage_factor = pixel_num / size(g,1);
        sensor_size = pixel_num * detector_size;
        deltaX = detector_size * shrinkage_factor;
        deltaY = detector_size * shrinkage_factor;

        if pad_size > 0
            g = padarray(g, [pad_size pad_size]);
        end
        [nx, ny] = size(g);

        Nx = nx;
        Ny = ny * nz * 2;   % 跟你原 TV 代码保持一致
        Nz = 1;

        E0 = ones(nx,ny);

        [Phase3D, Pupil] = MyMakingPhase3D(nx,ny,nz,lambda_um, ...
                                           deltaX,deltaY,deltaZ, ...
                                           offsetZ,sensor_size);

        E = MyFieldsPropagation(E0,nx,ny,nz,Phase3D,Pupil);

        g_vec = MyC2V(g(:));

        A  = @(f_twist) MyForwardOperatorPropagation(f_twist,E,nx,ny,nz,Phase3D,Pupil);
        AT = @(gmeas)   MyAdjointOperatorPropagation(gmeas,E,nx,ny,nz,Phase3D,Pupil);

        Psi = @(f,th) MyTVpsi(f,th,0.05,piter,Nx,Ny,Nz);
        Phi = @(f)    MyTVphi(f,Nx,Ny,Nz);

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

        f_c = MyV2C(f_vec);
        f_c = reshape(f_c, nx, ny, nz);

        mag = abs(f_c);

        % 去掉 pad
        if pad_size > 0
            mag = mag( pad_size+1:end-pad_size, ...
                       pad_size+1:end-pad_size, : );
        end

        maxVal = max(mag(:));
        if maxVal > 0
            vol = mag / maxVal;
        else
            vol = mag;
        end
    end

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

        [ny, nx] = size(imgNet);
        nz = numel(zVec);
        vol = repmat(sharp, [1, 1, nz]);
    end

    function [coords3D, axialCurves, roiBoxes, bw, img2D, stats] = localizeParticles3D(volData, zVecUm)
        img2D = max(volData, [], 3);
        img2D = normalizeImage(img2D);

        level = 0.7;
        bw = imbinarize(img2D, level);
        bw = bwareaopen(bw, 5);
        bw = imfill(bw, 'holes');

        stats = regionprops(bw, img2D, 'Centroid', 'WeightedCentroid', ...
                            'EquivDiameter', 'BoundingBox');
        if isempty(stats)
            coords3D = zeros(0, 3);
            axialCurves = {};
            roiBoxes = zeros(0, 4);
            return;
        end

        [nRows, nCols, nZ] = size(volData);
        coords3D = zeros(numel(stats), 3);
        axialCurves = cell(numel(stats), 1);
        roiBoxes = zeros(numel(stats), 4);

        for n = 1:numel(stats)
            centerXY = stats(n).WeightedCentroid;
            if any(~isfinite(centerXY))
                centerXY = stats(n).Centroid;
            end

            cx = round(centerXY(1));
            cy = round(centerXY(2));
            roiRadiusPx = max(3, ceil(stats(n).EquivDiameter / 2));

            x1 = max(1, cx - roiRadiusPx);
            x2 = min(nCols, cx + roiRadiusPx);
            y1 = max(1, cy - roiRadiusPx);
            y2 = min(nRows, cy + roiRadiusPx);
            roiBoxes(n, :) = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];

            axialCurve = zeros(nZ, 1);
            for k = 1:nZ
                roiSlice = volData(y1:y2, x1:x2, k);
                axialCurve(k) = mean(roiSlice(:));
            end

            [~, maxIdx] = max(axialCurve);
            focusSlice = volData(:, :, maxIdx);
            stats(n).EquivDiameter = estimateFocusDiameter( ...
                focusSlice, x1, x2, y1, y2, cx, cy, level, stats(n).EquivDiameter);
            coords3D(n, :) = [centerXY(1), centerXY(2), zVecUm(maxIdx)];
            axialCurves{n} = axialCurve;
        end
    end

    function focusDiamPx = estimateFocusDiameter(focusSlice, x1, x2, y1, y2, cx, cy, xyThreshold, fallbackDiamPx)
        roiFocus = focusSlice(y1:y2, x1:x2);
        roiFocus = normalizeImage(roiFocus);
        if ~any(roiFocus(:) > 0)
            focusDiamPx = fallbackDiamPx;
            return;
        end

        if exist('graythresh', 'file') == 2
            localLevel = graythresh(roiFocus);
        else
            localLevel = xyThreshold;
        end
        localLevel = max(0.35, min(0.85, max(localLevel, 0.6 * xyThreshold)));

        bwFocus = imbinarize(roiFocus, localLevel);
        bwFocus = bwareaopen(bwFocus, 1);
        bwFocus = imfill(bwFocus, 'holes');
        if ~any(bwFocus(:))
            focusDiamPx = fallbackDiamPx;
            return;
        end

        cc = bwconncomp(bwFocus, 8);
        localStats = regionprops(cc, roiFocus, 'EquivDiameter', 'WeightedCentroid', 'Centroid');

        cxLocal = min(size(bwFocus, 2), max(1, cx - x1 + 1));
        cyLocal = min(size(bwFocus, 1), max(1, cy - y1 + 1));
        centerIdx = sub2ind(size(bwFocus), cyLocal, cxLocal);

        pickIdx = [];
        for k = 1:cc.NumObjects
            if any(cc.PixelIdxList{k} == centerIdx)
                pickIdx = k;
                break;
            end
        end

        if isempty(pickIdx)
            centers = reshape([localStats.WeightedCentroid], 2, []).';
            if isempty(centers) || any(~isfinite(centers(:)))
                centers = reshape([localStats.Centroid], 2, []).';
            end
            dist2 = (centers(:, 1) - cxLocal).^2 + (centers(:, 2) - cyLocal).^2;
            [~, pickIdx] = min(dist2);
        end

        focusDiamPx = localStats(pickIdx).EquivDiameter;
        if ~(isfinite(focusDiamPx) && focusDiamPx > 0)
            focusDiamPx = fallbackDiamPx;
        end
    end

    function img = normalizeImage(img)
        img = double(img);
        imgMin = min(img(:));
        imgMax = max(img(:));
        if imgMax > imgMin
            img = (img - imgMin) / (imgMax - imgMin);
        else
            img = zeros(size(img));
        end
    end

    function hParticleFig = showParticleSpheres(coords3D, diams_um, pix_um)
        if isempty(coords3D)
            hParticleFig = [];
            return;
        end

        coordsUm = [coords3D(:,1) * pix_um, coords3D(:,2) * pix_um, coords3D(:,3)];
        if isempty(diams_um)
            baseRadius = max(40, 6 * pix_um);
        else
            baseRadius = max([40, 6 * pix_um, median(diams_um) * 1.5]);
        end

        hParticleFig = figure( ...
            'Name', 'Sparse Particle Field', ...
            'NumberTitle', 'off', ...
            'Color', [0.05 0.07 0.12], ...
            'Position', [1150 120 820 640]);
        ax = axes('Parent', hParticleFig);
        hold(ax, 'on');

        [sx, sy, sz] = sphere(20);
        colors = lines(max(size(coordsUm, 1), 7));
        for n = 1:size(coordsUm, 1)
            radiusUm = max(baseRadius, diams_um(min(n, numel(diams_um))) * 1.2 / 2);
            surf(ax, ...
                sx * radiusUm + coordsUm(n, 1), ...
                sy * radiusUm + coordsUm(n, 2), ...
                sz * radiusUm + coordsUm(n, 3), ...
                'FaceColor', colors(mod(n-1, size(colors,1)) + 1, :), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.95);
        end

        grid(ax, 'on');
        axis(ax, 'equal');
        view(ax, 32, 24);
        xlabel(ax, ['x (' char(181) 'm)']);
        ylabel(ax, ['y (' char(181) 'm)']);
        zlabel(ax, ['z (' char(181) 'm)']);
        title(ax, 'Localized Particle Field');
        set(ax, 'Color', [0.10 0.11 0.15], ...
            'XColor', [0.85 0.90 1.00], ...
            'YColor', [0.85 0.90 1.00], ...
            'ZColor', [0.85 0.90 1.00], ...
            'FontName', 'Times New Roman');
        camlight(ax, 'headlight');
        lighting(ax, 'gouraud');
        drawnow;
    end
% =========================================================================
    % [新增] 简易粒径统计回调 (复用主界面数据)
    % =========================================================================
    function onSimpleStats(~, ~)
        handles = guidata(hFig); % 获取主界面所有数据
        
        % 1. 数据校验：直接使用主界面共享的 VolumeData
        if isempty(handles.VolumeData)
            errordlg('请先在主界面点击"计算并绘制3D图"生成数据！', '无数据');
            return;
        end
        
        % 2. 参数获取：直接使用主界面共享的控件参数
        pix_str = get(handles.hPixel, 'String');
        pix_um = str2double(pix_str);
        if isnan(pix_um), pix_um = 1; end % 默认值防错
        
        % 3. 图像处理 (后台自动完成)
        set(handles.hStatusText, 'String', 'Status: Calculating particle stats...'); drawnow;
        
        % MIP 投影：将 3D 数据压平为 2D
        volData = handles.VolumeData;
        img2D = max(volData, [], 3);
        
        % 归一化 & 自动二值化
        img2D = (img2D - min(img2D(:))) / (max(img2D(:)) - min(img2D(:)));
        %level = graythresh(img2D);
       
        
        level=0.7;
       
        
        
        
        bw = imbinarize(img2D, level);
        bw = bwareaopen(bw, 5);   % 滤除 < 5像素的噪点
        bw = imfill(bw, 'holes'); % 填充内部空洞
        
        % 4. 统计计算
        stats = regionprops(bw, 'Centroid', 'EquivDiameter');
        
        if isempty(stats)
            msgbox('当前阈值下未检测到颗粒，请检查图像质量。', '统计结果');
            set(handles.hStatusText, 'String', 'Status: No particles found');
            return;
        end
        
        diams_px = [stats.EquivDiameter];
        diams_um = diams_px * pix_um; % 转换为微米
        
        % 5. 结果展示 (只弹出一个报告窗口)
        hReport = figure('Name', '粒径统计结果报告', ...
            'NumberTitle', 'off', ...
            'Position', [300 300 800 400], ... % 宽一点，并排显示
            'Color', 'white');
            
        % 左图：显示哪里被识别成了颗粒 (画红圈)
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
        ylabel(ax2, '颗粒数量');
        
        % 底部：文字汇总
        infoStr = sprintf(['统计汇总:  ', ...
            '颗粒总数: %d 个   |   ', ...
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
            
        set(handles.hStatusText, 'String', 'Status: Statistics report generated');
    end

    function onSimpleStatsV2(~, ~)
        handles = guidata(hFig);

        if isempty(handles.VolumeData)
            errordlg('请先在主界面完成 3D 重建。', '无数据');
            return;
        end
        if isempty(handles.VolumeZVecUm)
            errordlg('当前结果缺少 z 轴采样信息，请重新重建。', '无数据');
            return;
        end

        pix_um = str2double(get(handles.hPixel, 'String'));
        if isnan(pix_um)
            pix_um = 1;
        end

        set(handles.hStatusText, 'String', 'Status: Localizing particles in 3D...'); drawnow;

        volData = handles.VolumeData;
        zVecUm = handles.VolumeZVecUm(:);
        [coords3D, axialCurves, roiBoxes, ~, img2D, stats] = localizeParticles3D(volData, zVecUm);

        if isempty(stats)
            msgbox('当前阈值下未检测到粒子，请检查重建结果。', '统计结果');
            set(handles.hStatusText, 'String', 'Status: No particles found');
            return;
        end

        diams_px = [stats.EquivDiameter];
        diams_um = diams_px * pix_um;
        zStepUm = 1;
        if numel(zVecUm) > 1
            zStepUm = median(abs(diff(zVecUm)));
        end
        if ~(isfinite(zStepUm) && zStepUm > 0)
            zStepUm = 1;
        end
        zSpanUm = max(zStepUm, abs(zVecUm(end) - zVecUm(1)) + zStepUm);
        reconVolumeUm3 = size(volData, 1) * size(volData, 2) * (pix_um ^ 2) * zSpanUm;
        numConcPerMm3 = size(coords3D, 1) * 1e12 / max(reconVolumeUm3, eps);

        diamVecUm = diams_um(:);
        volWeights = diamVecUm .^ 3;
        if any(volWeights > 0)
            [diamSorted, sortIdx] = sort(diamVecUm);
            cumFrac = cumsum(volWeights(sortIdx)) / sum(volWeights);
            idx50 = find(cumFrac >= 0.5, 1, 'first');
            if isempty(idx50)
                mvdUm = diamSorted(end);
            elseif idx50 == 1
                mvdUm = diamSorted(1);
            else
                prevFrac = cumFrac(idx50 - 1);
                nextFrac = cumFrac(idx50);
                prevDiam = diamSorted(idx50 - 1);
                nextDiam = diamSorted(idx50);
                if nextFrac > prevFrac
                    interpT = (0.5 - prevFrac) / (nextFrac - prevFrac);
                else
                    interpT = 0;
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
        handles.ParticleCoords3D = coords3D;
        handles.ParticleAxialCurves = axialCurves;
        guidata(hFig, handles);

        hReport = figure('Name', '粒子三维定位报告', ...
            'NumberTitle', 'off', ...
            'Position', [220 180 1100 520], ...
            'Color', 'black');

        ax1 = subplot(1, 3, 1);
        imshow(img2D, []); hold(ax1, 'on');
        set(ax1, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
        title(ax1, 'MIP 与轴向搜索 ROI', 'Color', 'w');
        centers = coords3D(:, 1:2);
        radii = diams_px' / 2;
        viscircles(ax1, centers, radii, 'Color', 'r', 'LineWidth', 0.5);
        for n = 1:size(roiBoxes, 1)
            rectangle(ax1, 'Position', roiBoxes(n, :), 'EdgeColor', 'g', 'LineWidth', 0.8);
            text(ax1, centers(n, 1), centers(n, 2), sprintf('%d', n), ...
                'Color', 'y', 'FontSize', 8, 'HorizontalAlignment', 'center');
        end

        ax2 = subplot(1, 3, 2);
        histogram(ax2, diams_um, 20, 'FaceColor', [0.2 0.6 0.8]);
        set(ax2, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
        grid(ax2, 'on');
        title(ax2, '粒径分布', 'Color', 'w');
        xlabel(ax2, ['直径 (' char(181) 'm)'], 'Color', 'w');
        ylabel(ax2, '颗粒数量', 'Color', 'w');

        ax3 = subplot(1, 3, 3);
        hold(ax3, 'on');
        set(ax3, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
        curvePeaks = zeros(numel(axialCurves), 1);
        for n = 1:numel(axialCurves)
            curvePeaks(n) = max(axialCurves{n});
        end
        [~, order] = sort(curvePeaks, 'descend');
        for n = 1:min(5, numel(order))
            idxCurve = order(n);
            plot(ax3, zVecUm, axialCurves{idxCurve}, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('P%d: z=%.0f %sm', idxCurve, coords3D(idxCurve, 3), char(181)));
        end
        grid(ax3, 'on');
        xlabel(ax3, ['z (' char(181) 'm)'], 'Color', 'w');
        ylabel(ax3, 'Mean ROI Intensity', 'Color', 'w');
        title(ax3, '轴向光强曲线 I(z)', 'Color', 'w');
        hLgd = legend(ax3, 'Location', 'best');
        set(hLgd, 'TextColor', 'w', 'Color', 'k', 'EdgeColor', [0.6 0.6 0.6]);

        infoStr = sprintf(['统计汇总  粒子总数: %d 个  |  平均直径: %.2f %sm  |  ', ...
            'z 范围: %.0f~%.0f %sm  |  步长: %.0f %sm'], ...
            length(diams_um), mean(diams_um), char(181), ...
            zVecUm(1), zVecUm(end), char(181), median(diff(zVecUm)), char(181));

        uicontrol('Parent', hReport, 'Style', 'text', ...
            'String', infoStr, ...
            'Units', 'normalized', ...
            'Position', [0 0.92 1 0.08], ...
            'BackgroundColor', 'black', ...
            'ForegroundColor', 'white', ...
            'FontSize', 10, 'FontWeight', 'bold');

        metricStr = sprintf(['中值体积直径: %.2f %sm  |  平均有效直径 D[3,2]: %.2f %sm  |  ', ...
            '粒子数浓度: %.3g 个/mL'], ...
            mvdUm, char(181), meanEffDiamUm, char(181), numConcPerMm3);
        uicontrol('Parent', hReport, 'Style', 'text', ...
            'String', metricStr, ...
            'Units', 'normalized', ...
            'Position', [0 0.875 1 0.04], ...
            'BackgroundColor', 'black', ...
            'ForegroundColor', [0.82 0.90 1.00], ...
            'FontSize', 9, ...
            'HorizontalAlignment', 'left');

        coordData = [(1:size(coords3D, 1))', coords3D, diams_um(:)];
        uitable('Parent', hReport, ...
            'Units', 'normalized', ...
            'Position', [0.02 0.01 0.96 0.18], ...
            'Data', coordData, ...
            'ColumnName', {'ID', 'x (px)', 'y (px)', ['z (' char(181) 'm)'], ['d (' char(181) 'm)']}, ...
            'BackgroundColor', [0 0 0], ...
            'ForegroundColor', [1 1 1]);

        hParticleFig = showParticleSpheres(coords3D, diams_um, pix_um);
        if ~isempty(hParticleFig) && isgraphics(hParticleFig)
            figure(hParticleFig);
        end

        set(handles.hStatusText, 'String', sprintf('Status: %d particles localized in 3D', size(coords3D, 1)));
    end

    function onRedrawV2(~, ~)
        handles = guidata(hFig);

        set(handles.hStatusText, 'String', 'Status: Redrawing with current alpha...');
        drawnow;

        if isempty(handles.VolumeData)
            errordlg('请先点击“计算并绘制3D图”生成体数据', '提示');
            set(handles.hStatusText, 'String', 'Status: Redraw failed (no volume data)');
            drawnow;
            return;
        end

        onRedraw([], []);
        set(handles.hStatusText, 'String', 'Status: Redraw finished');
        drawnow;
    end
end
