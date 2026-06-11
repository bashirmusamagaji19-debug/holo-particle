%DINGWEI2(1) 定位原型脚本——基于最大强度投影和阈值分割的3D粒子定位
%
%   【用途】
%   从3D重建体数据中定位粒子位置。核心流程：
%   1. 计算最大强度投影（MIP）；
%   2. 阈值过滤弱信号（10%最大值以下置零）；
%   3. 二值化+连通域分析筛选有效粒子；
%   4. 在原始3D数据中找到每个粒子的z坐标（最大强度层）；
%   5. 生成HSV彩色编码图（H=深度，V=强度）；
%   6. 显示各层重建切片。
%
%   【在工程中的位置】
%   true_way目录下的历史定位原型脚本，不在主链路中。
%   这是最早期的粒子定位方法：基于阈值+连通域的简单方案。
%   当前主GUI使用更完善的定位算法（如3D峰值检测+亚像素拟合）。
%
%   【注意】
%   文件名含括号和数字，为版本迭代存档。代码中部分中文注释为乱码
%   （原始编码问题），已尽量还原其含义。

clear all;
close all;
clc;

%% ---- 加载重建数据 ----
% load insect9-2-deltaZ200-result-bp.mat      % 早期数据集
load f_reconstruct_3.mat     % 重建结果体数据
load data.mat                % 原始数据

% 显示重建体数据的3D视图（permute调整轴顺序）
figure,show3d(permute(abs(f_reconstruct),[2,3,1]),0.7)

% f= transf; % 
f = abs(data);  % 取数据幅度作为工作体数据

%% ---- 获取数据尺寸 ----
[Nx,Ny,Nz] = size(f);

%% ---- 最大强度投影 + 弱信号滤除 ----
o = zeros(Nx,Ny,Nz);        % 预分配输出体数据
a = zeros(Nx,Ny);            % 预分配MIP图
for i = 1:Nx
    for j = 1:Ny
        a(i,j) = max(abs(f(i,j,:)));   % 沿z方向取最大值→MIP
        for k = 1:Nz
            if abs(f(i,j,k)) <= a(i,j)*0.1   % 强度低于该像素MIP的10%则置零
                f(i,j,k) = 0;
            end   
        end  
    end
end

% 显示最大强度投影
figure()
imshow(a,[]);
title('最大强度投影')

%% ---- 二值化与连通域分析 ----
% 阈值：取MIP最大值的85%作为二值化阈值
ax1 = a > 0.85*max(max(a));

% 对二值图进行连通域标记
imLabel = bwlabel(ax1);                     % 各连通域标记
stats = regionprops(imLabel,'Area');         % 计算各连通域面积
area = cat(1,stats.Area);
% 只保留面积>=3像素的连通域（去除孤立噪点）
imLabel = ismember(imLabel, find([stats.Area] >= 3));

figure(),imshow(imLabel,[]);
title('二值化图')

%% ---- 在MIP上提取筛选后的粒子区域 ----
ax2 = imLabel .* a;  % 将二值掩模与MIP强度相乘，保留粒子区域的强度值
figure()
imshow(ax2,[]);
title('选出的粒子区域')

%% ---- 保存ground truth图像 ----
load counter_img.mat
savename = sprintf('./dataset_CS/groundtruth/gt%g.jpeg',counter_img);

%% ---- 在3D体数据中定位每个粒子的z坐标 ----
h1 = zeros(Nx,Ny);  % 预分配z坐标图（深度图）
for i = 1:Nx
    for j = 1:Ny
        if ax2(i,j) == 0
            ax2(i,j) = 0;
        else
            % 在第(i,j)像素的z方向上，找到强度等于MIP值的层索引
            e = find(ax2(i,j) == abs(f(i,j,:)));   % 该像素的z位置
            o(i,j,e) = ax2(i,j);                    % 在3D体数据中标记粒子位置
            h1(i,j) = e;                             % 记录z坐标
        end
    end
end

%% ---- 计算平均z深度 ----
[x,y] = find(ax1~=0);  % 找到所有粒子像素的坐标
t = size(x,1);
sum_val = 0;
for i = 1:t
    sum_val = sum_val + h1(x(i),y(i));
end
avg_z = sum_val/t   % 显示平均z深度

%% ---- 生成HSV彩色编码图 ----
% H（色调）= z深度归一化值，S（饱和度）= 1，V（明度）= 强度
% 这样不同深度的粒子显示不同颜色
figure,imshow(h1,[]);
colorbar
N = size(o,1);
s = ones(N);
hsv = ones(N,N,3);
hsv(:,:,1) = h1 / max(max(ax2));  % H通道：归一化深度
hsv(:,:,2) = s;                    % S通道：满饱和
hsv(:,:,3) = ax2;                  % V通道：强度值
k = hsv2rgb(hsv);                  % 转换为RGB
figure,imshow(k,[]);
colorbar

%% ---- 显示各层切片 ----
figure;
counter = 1;
for i = 1:size(o,3)
    subplot(6,6,counter)
    imshow(o(:,:,i),[])
    counter = counter+1;
    title(num2str(i))
    if counter>36
        figure;
        counter = 1;
    end
end
