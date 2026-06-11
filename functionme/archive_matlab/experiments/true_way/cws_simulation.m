function [hlo,xt]=cws_simulation(z1,obj,pix_pitch,lambda)
%CWS_SIMULATION 卷积传播仿真函数
%
%   【用途】
%   基于卷积法（Fresnel卷积）实现全息图的正向/反向传播仿真。
%   给定物体复振幅和传播距离，计算传播后的复振幅和模拟全息图。
%   包含多个辅助子函数，覆盖Fresnel方程法、卷积法、角谱法等传播方式。
%
%   【在工程中的位置】
%   true_way目录下的历史原型函数。被BP_G.m调试脚本调用。
%   主链路不使用此函数（主链路走AS_reconstructAngularSpectrum系列）。
%   本函数是早期卷积传播的实验实现，用于验证传播算法的正确性。
%
%   【输入参数】
%       z1        : 传播距离（米），负值表示反向传播（重建）
%       obj       : 输入图像/物体振幅（2D矩阵）
%       pix_pitch : 像素间距（米）
%       lambda    : 波长（米）
%
%   【输出参数】
%       hlo : 传播后的复振幅
%       xt  : 原始物体（=obj的副本）

%% ==================== 注释掉的bin文件读取代码 ====================
% 以下大段注释为早期从bin格式实验数据中提取全息图和参考光的代码，
% 保留作为参考，当前不再使用。

%% extract image from bin files start
% % % % % % % % % % num_sets = 120; 
% % % % % % % % % % num_sets_ref=10;
% % % % % % % % % % num_per_set=3; 
% % % % % % % % % % height=1280; 
% % % % % % % % % % width=1024; 
% % % % % % % % % % data = zeros(height,width,num_per_set,num_sets); 
% % % % % % % % % % for i=1:num_sets 
% % % % % % % % % % % % % %     fp=fopen(sprintf('.\\g1_PST_exp3_201712111522\\img_%03d.bin',i),'rb');
% % % % % % % % % %     fp=fopen(sprintf('.\\PST_exp3_201801241642\\img_%03d.bin',i),'rb'); 
% % % % % % % % % %     data(:,:,:,i)=reshape(fread(fp,height*width*num_per_set,'uint8'),height,width,num_per_set); 
% % % % % % % % % %     fclose(fp); 
% % % % % % % % % % end
%% extract image from bin files end
%% save the second image of each bin file start
% % % % % % % % % % for i=1:num_sets
% % % % % % % % % %     hologr=rot90(data(:,:,2,i),-1);
% % % % % % % % % % % % % %     filename=sprintf('.\\g1_PST_exp3_201712111522\\pholograms\\hologram_%d.bmp',i);
% % % % % % % % % %     filename=sprintf('.\\PST_exp3_201801241642\\pholograms\\hologram_%d.bmp',i);
% % % % % % % % % %     imwrite(uint8(hologr),filename);
% % % % % % % % % % end
%% save the second image of each bin file end

%% extract the reference start
% % % % % % % % % % % for i=1:num_sets_ref 
% % % % % % % % % % %     fp=fopen(sprintf('.\\PST_exp3_201801241647\\img_%03d.bin',i),'rb'); 
% % % % % % % % % % %     data_ref(:,:,:,i)=reshape(fread(fp,height*width*num_per_set,'uint8'),height,width,num_per_set); 
% % % % % % % % % % %     fclose(fp); 
% % % % % % % % % % % end
%% extract the reference end

%% save the second image of each bin file start
% % % % % % % % % % for i=1:num_sets_ref
% % % % % % % % % %     hologr=rot90(data_ref(:,:,2,i),-1);
% % % % % % % % % %     filename=sprintf('.\\PST_exp3_201801241647\\reference\\reference_%d.bmp',i);
% % % % % % % % % %     imwrite(uint8(hologr),filename);
% % % % % % % % % % end
%% save the second image of each bin file end
% % % % height=1280; 
% % % % width=1024; 

%% ==================== 主传播计算 ====================
% obj1=double(imread('brain','BMP')); % recorded hologram
xt = obj;                             % 保存原始物体副本
[objx, objy] = size(obj);            % 获取物体尺寸

%parameters display (start)
um = 1e-6;
mm = 1e-3;                
k = 2*pi/lambda;                     % 波数 k = 2π/λ

%% ---- 计算物体物理尺寸 ----
obj_dime_x = objx * pix_pitch;       % x方向物理尺寸（米）
obj_dime_y = objy * pix_pitch;       % y方向物理尺寸（米）

%% ---- 卷积法频域采样 ----
% 生成频域坐标网格，用于卷积法传播
[f_x, f_y] = sampl_freq(objx, objy, obj_dime_x, obj_dime_y);

%% ---- 执行Fresnel卷积传播 ----
[hlo] = FrPr_conlu(obj, lambda, z1, f_x, f_y);
hlo = conj(hlo);                     % 取共轭（反向传播需要）

%% ---- 模拟全息图记录 ----
% g = 1 + |hlo|² + hlo* + hlo  模拟离轴全息图的四项
g = 1 + abs(hlo).^2 + conj(hlo) + hlo;  % 模拟相机记录的全息图
% % % % figure; imagesc(hlo_abs); colormap gray; truesize; title('the amplitude of the hologram');
% % % % gg=double(imread('capt_hologram','BMP')); % recorded hologram
% % % % 
% % % % [re_im]=FrPr_conlu(gg-1, lambda, -z1, f_x, f_y);
% % % % re_im=abs(re_im)*255/max(max(abs(re_im)));
% % % % 
% % % % figure; imagesc(re_im); colormap gray; truesize; title('reconstructed image');

%% ---- 归一化并转为uint8 ----
g_temp = abs(g).^2;                  % 取强度
g_temp = g_temp * 255 / max(max(g_temp));  % 归一化到0-255
cap_holo = uint8(g_temp);            % 转为8位无符号整数
% imwrite(cap_holo,'capt_hologram.bmp');
% % % % % [PSNR] = psnr(re_im, obj);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ==================== 辅助子函数 ====================

function [OBJ, X,Y] = pic_read(name, M, N)
%PIC_READ 读取图片并居中嵌入到M×N画布
%   M,N为目标画布尺寸；若M<=0或N<=0则使用原始尺寸
%   X,Y为原始图像尺寸
   OBJ_org = imread(name);
   OBJ_info = imfinfo(name);
   if(OBJ_info.BitDepth>8)
       OBJ_org = double(rgb2gray(OBJ_org));  % 高位深先转灰度
   else
        OBJ_org = double(OBJ_org);
   end
     [Y,X] = size(OBJ_org);
   
   if M>0 &&N>0
       % 居中嵌入到M×N画布
       OBJ = zeros(M, N);
       OBJ( N/2-Y/2+1 : N/2+Y/2, N/2-X/2+1 : N/2+X/2) = OBJ_org;
   else
       OBJ = OBJ_org;  % 保持原始尺寸
   end
   
    figure();
    imshow(OBJ,[]);
end

function [sa_x,sa_y]=sampling(in_x, in_y, X, Y)
%SAMPLING 方程法采样网格生成
%   生成空域坐标网格，用于Fresnel方程法传播
[sa_x,sa_y]=meshgrid(-X/2*in_x : in_x : -X/2*in_x+(X-1)*in_x, -Y/2*in_y : in_y : -Y/2*in_y+(Y-1)*in_y);
end

function [hlo] = FrPr(obj,sigma,sa_x_obj,sa_y_obj,sa_x_hlo, sa_y_hlo)
%FrPr Fresnel方程法正向传播（含fftshift）
%   通过二次相位因子和FFT实现Fresnel传播
u1 = obj.*exp(1i*sigma*(sa_x_obj.^2+sa_y_obj.^2));  % 乘以物面二次相位
u1 = fftshift(fft2(u1));                              % FFT并中心化
hlo = u1.*exp(1i*sigma*(sa_x_hlo.^2+sa_y_hlo.^2));   % 乘以全息面二次相位
end

function [hlo] = FrPr1(obj,sigma,sa_x_obj,sa_y_obj,sa_x_hlo, sa_y_hlo)
%FrPr1 Fresnel方程法正向传播（无fftshift版本）
u1 = obj.*exp(1i*sigma*(sa_x_obj.^2+sa_y_obj.^2));
u1 = (fft2(u1));                                      % FFT不做shift
hlo = u1.*exp(1i*sigma*(sa_x_hlo.^2+sa_y_hlo.^2));
end

function [re_im] = In_FrPr(hlo,sigma,sa_x_hlo,sa_y_hlo,sa_x_objre,sa_y_objre)
%In_FrPr Fresnel方程法反向传播（重建）
u2 = hlo.*exp(-1i*sigma*(sa_x_hlo.^2+sa_y_hlo.^2));  % 乘以全息面逆二次相位
u2 = ((ifft2(u2)));                                    % IFFT
re_im = u2.*exp(-1i*sigma*(sa_x_objre.^2+sa_y_objre.^2)); % 乘以物面逆二次相位
end





function []=im_s(hlo) 
%IM_S 显示复振幅的幅度图像
im_temp = abs(hlo);
imshow(im_temp/max(max(im_temp)),[]);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 卷积法频域采样
function [f_x, f_y]=sampl_freq(X, Y, x_dimen, y_dimen)
%SAMPL_FREQ 生成卷积法所需的频域采样网格
%   X,Y为图像像素尺寸；x_dimen,y_dimen为物理尺寸
%   频率间隔 = 1/物理尺寸
[f_x,f_y]=meshgrid((-X/2)*(1/x_dimen):(1/x_dimen):(X/2-1)*(1/x_dimen),(-Y/2)*(1/y_dimen):(1/y_dimen):(Y/2-1)*(1/y_dimen));
end

%% Fresnel卷积法传播（主用函数）
function [hlo]=FrPr_conlu(obj, lambda, d, f_x, f_y)
%FrPr_CONLU Fresnel卷积法传播
%   lambda: 波长
%   obj: 物体复振幅
%   d: 传播距离
%   f_x, f_y: 频域坐标网格（由sampl_freq生成）
%
%   传递函数：H = exp(i·k·d)·exp(-i·π·λ·d·(fx²+fy²))
%   这是Fresnel近似下的卷积传播核

% 计算频率间隔，用于确定传播通带半径
delta_fx = abs(f_x(1,1)-f_x(1,2));
R = ((1/(lambda*d))+delta_fx.^2)/(2*delta_fx);  % 通带半径
mask = (f_x.^2+f_y.^2<R.^2);                     % 通带掩模（当前未使用）

pi*lambda*d;
k = 2*pi/lambda;

% Fresnel卷积传递函数（无通带掩模版本）
H = exp(1i*k*d)*exp(-1i*pi*lambda*d*(f_x.^2+f_y.^2)); % without mask
% H = exp(1i*k*d)*exp(-1i*pi*lambda*d*(f_x.^2+f_y.^2)).*mask; % with mask（带通带限制）

H_temp = ifftshift(H);             % 将H中心移到角落以匹配FFT输出
O = fft2(ifftshift(obj));          % 物体FFT（先ifftshift使DC到角落）
hlo_temp = H_temp .* O;            % 频域相乘
hlo = fftshift(ifft2(hlo_temp));   % IFFT后fftshift恢复中心
end

%% 角谱法核函数生成（备用，当前未在主流程中调用）
function [kernals]=kernal_making(lambda, d, f_x, f_y)
%KERNAL_MAKING 生成角谱法传递函数核
%   H = exp(i·k·d·sqrt(1 - λ²fx² - λ²fy²))
%   这是严格角谱传播（非Fresnel近似）
k = 2*pi/lambda;
sqrt_temp = sqrt(1-(lambda*f_x).^2-(lambda*f_y).^2);
H = exp(1i*k*d.*sqrt_temp);
kernals = fftshift(H);
end

%% 角谱法反向传播（备用）
function [hlo]=FrPr_incon_pt(obj, lambda, d, f_x, f_y)
%FrPr_INCON_PT 角谱法传播（严格版，非Fresnel近似）
%   使用精确的角谱传递函数进行传播
k = 2*pi/lambda;

% 严格角谱传递函数
sqrt_temp = sqrt(1-(lambda*f_x).^2-(lambda*f_y).^2);
H = exp(1i*k*d.*sqrt_temp);

H_temp = fftshift(H);
O = fft2(fftshift(obj));
hlo_temp = H_temp .* O;
hlo = ifftshift(ifft2(hlo_temp));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% PSNR计算
function [PSNR] = psnr(f1, f2)    
%PSNR 计算峰值信噪比
%   用于评估重建图像与原始图像之间的质量
f1 = im2uint8(f1/max(max(f1)));
f2 = im2uint8(f2/max(max(f2)));
k = 8;       % 图像位深
fmax = 2.^k - 1; 
a = fmax.^2; 
e = double(f1) - double(f2); 
[m, n] = size(e); 
b = sum(sum(e.^2)); 
PSNR = 10*log10(m*n*a/b);
disp(PSNR);
end

%% 相对误差计算
function [re_error]=relative_error(re,or)
%RELATIVE_ERROR 计算相对误差
re_error = sum(sum((re-or).^2))/sum(sum(or.^2));
end

%% 相位归一化到[0,2π]
function [B]=angle_2pi(A)
%ANGLE_2PI 将相位值归一化到[0, 2π]区间
[Y,X] = size(A);
exp_A = exp(1i*A);
angle_exp_A = angle(exp_A);
I = find(angle_exp_A<=0);
angle_exp_A(I) = 2*pi - angle_exp_A(I)*(-1);  % 负相位翻转为正值
B = reshape(angle_exp_A,X,Y);
end

%% 随机置乱图像
function [obj_n,Inx]=en_rand_image(OBJ,enr)
%EN_RAND_IMAGE 按索引enr对图像行进行随机置乱
[y,x] = size(OBJ);
en_im = reshape(OBJ(enr),x,y);
obj_n = en_im;
[ori,Inx] = sort(enr);  % Inx记录排序索引，用于后续还原
end

%% 随机置乱还原图像
function [ori]=de_rand_image(re_im,Inx)
%DE_RAND_IMAGE 按索引Inx还原被置乱的图像
[y,x] = size(re_im);
ori = re_im(Inx);
ori = reshape(ori,x,y);
end

%% 迭代过程图像处理（用于相位恢复等迭代算法）
function [real_n,imag_n]=iteritative_process_image(OBJ,num)
%ITERITATIVE_PROCESS_IMAGE 对图像进行迭代随机相位处理
%   将实部和虚部分离，每步叠加随机相位
[y,x] = size(OBJ);
real_n = OBJ;
for n = 1:num
    alpha_exp = exp(1i*2*pi*rand(x,y));  % 随机相位因子
    obj_ph = real_n .* alpha_exp;         % 乘以随机相位
     if n>1
        imag_n = imag_n .* alpha_exp;
    end
    A = rand(x,y);                        % 随机阈值
    obj_ph_real = real(obj_ph) - A;       % 实部减去随机值
    real_n = obj_ph_real;
    
    obj_ph_imag = imag(obj_ph)*1i + A;    % 虚部加随机值
   
     if n==1
         imag_n = obj_ph_imag;
     else
         imag_n = obj_ph_imag + imag_n;   % 累加虚部
    end
end
end

%% RMS误差计算
function [rms_e]=rms_error(re,or)
%RMS_ERROR 计算均方根误差
[Y,X] = size(re);
ms_e = sum(sum((re-or).^2))*(1/(X*Y));
rms_e = sqrt(ms_e);
end



% % % figure(); 
% % % imshow(temp_fo,[]);
% % % white_blob=~temp_fo;
% % % white_blob=imfill(white_blob,'holes');
% % % stats_final=regionprops('table', white_blob, 'Centroid', 'MajorAxisLength', 'MinorAxisLength','EquivDiameter','Area','Image');
% % % centers = stats_final.Centroid;
% % % diameters = mean([stats_final.MajorAxisLength stats_final.MinorAxisLength],2);
% % % radii = diameters/2;
% % % hold on
% % % viscircles(centers,(radii+8));
% % % hold off
