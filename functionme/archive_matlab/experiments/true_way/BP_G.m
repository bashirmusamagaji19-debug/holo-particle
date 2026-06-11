%BP_G 实验/调试脚本——全息图读取与多深度卷积传播重建
%
%   【用途】
%   从BMP文件读取彩色全息图，选取RGB单通道和感兴趣区域（ROI），
%   使用卷积传播法（cws_simulation）在不同深度进行扫描重建，
%   最终用show3d进行3D体绘制显示。
%
%   【在工程中的位置】
%   true_way目录下的历史调试脚本，不在主链路中。
%   这是最早期的原型验证流程：手动选区 → 逐层重建 → 3D显示。
%   主GUI已将这些步骤整合为自动化流程。
%
%   【注意】
%   硬编码了原始开发环境的文件路径，需根据实际情况修改。
%   mm = 1e-3 为毫米到米的换算因子。

clc
clear
close all
addpath("./function/")   % 添加函数目录到搜索路径

%% ---- 参数定义 ----
mm = 1e-3;                % 毫米→米换算因子
pix_pitch = 2.2e-6;       % 像素间距 2.2μm
lambda1 = 638e-9;         % 红光波长 638nm

%% ---- 读取全息图 ----
% 读取彩色全息图（路径为原始开发环境，需修改）
hlo = imread("C:\Users\HAO ZHITAO\MVS\Data\2023-5-12\Image_6.bmp");

% 分离RGB三通道
t1 = hlo(:,:,1);  % R通道
t2 = hlo(:,:,2);  % G通道
t3 = hlo(:,:,3);  % B通道

%% ---- 显示原始全息图及三通道 ----
figure()
subplot(2,2,1),imshow(hlo,[]);title("hologram")       % 原始彩色全息图
subplot(2,2,2),imshow(t1,[]);title("R-channel")        % R通道
subplot(2,2,3),imshow(t2,[]);title("G-channel")        % G通道
subplot(2,2,4),imshow(t3,[]);title("B-channel")        % B通道
t=1

%% ---- 手动选取感兴趣区域（ROI） ----
figure,imshow(t1,[])
h = imrect;  % 交互式矩形选区工具
% 拖动鼠标获得兴趣区域，pos有四个值：兴趣区域左上角的像素坐标和区域的长宽
position = getPosition(h);

% 从position中提取ROI参数
start_y = floor(position(2));       % ROI起始y坐标
start_x = floor(position(1));       % ROI起始x坐标
rect_size_y = floor(position(3));   % ROI高度
rect_size_x = floor(position(4));   % ROI宽度
max_rect = 1023;                    % 固定裁剪为1023×1023像素
rect_size_y = max_rect;
rect_size_x = max_rect;

% 显示裁剪后的彩色ROI（加40亮度偏移以增强显示）
figure,imshow(hlo(start_y:start_y+rect_size_y,start_x:start_x+rect_size_x,:)+40,[])

% 提取R通道的ROI区域用于后续重建
g = t1(start_y:start_y+rect_size_y,start_x:start_x+rect_size_x);
figure,imshow(t1(start_y:start_y+rect_size_y,start_x:start_x+rect_size_x),[])
figure,imshow(t2(start_y:start_y+rect_size_y,start_x:start_x+rect_size_x),[])

%% ---- 多深度卷积传播重建 ----
lambda2 = 520e-9;              % G通道波长520nm（此脚本中未实际使用lambda2）
dis_start = 2*mm;              % 起始深度 2mm
dis_interval = 0.010*mm;       % 深度步进 10μm
dis_end = 4.6*mm;              % 终止深度 4.6mm
figure;
pause on
cnt = 1;
% 逐深度扫描重建
for dis = dis_start:dis_interval:dis_end
    % 调用cws_simulation进行卷积传播，负号表示反向传播
    [img, obj] = cws_simulation(-1*dis, g, pix_pitch, lambda1);
    tt = max(max(abs(img))) - abs(img);  % 反转强度（使粒子呈亮点）
    imshow(tt,[])
    title([num2str(dis*1e3),'mm'])
    volume(:,:,cnt) = tt;  % 将当前层存入3D体数据
    cnt = cnt+1;
    pause(0.5)  % 每层暂停0.5秒，便于观察重建过程
end

%% ---- 3D体绘制显示 ----
figure,show3d(volume, 0.36, size(volume,1), size(volume,3));
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1.2);
xlabel('\fontname{Times New Roman}\fontsize{17}\it{x}/\fontname{宋体}\fontsize{17}\rm像素')
ylabel('\fontname{Times New Roman}\fontsize{17}\it{y}/\fontname{宋体}\fontsize{17}\rm像素')
zlabel('\fontname{Times New Roman}\fontsize{17}\it{z}/\fontname{宋体}\fontsize{17}\rm像素')
set(gca,'Position',[0.25 0.25 0.5 0.5]);
