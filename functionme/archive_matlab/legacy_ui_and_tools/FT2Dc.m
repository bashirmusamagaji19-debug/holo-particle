%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 二维中心化傅里叶变换 (2D Centered Fourier Transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 【用途】
%   对二维输入矩阵执行中心化傅里叶变换，使零频(DC)分量位于频谱中心，
%   而非默认的四个角。在全息重建和角谱传播中，中心化频谱便于直观分析
%   和滤波操作。
%
% 【数学公式】
%   标准FFT的零频位于左上角。中心化通过移位性质实现：
%     G(u,v) = (-1)^(m+n) · FFT{ (-1)^(m+n) · g(m,n) }
%   其中 (-1)^(m+n) = exp(jπ(m+n))，对输入和输出各乘一次该因子，
%   等效于 fftshift(fft2(fftshift(x)))，但计算效率更高（避免数据重排）。
%
% 【输入】
%   in  - 二维矩阵 (Nx × Ny)，可为实数或复数
%
% 【输出】
%   out - 二维矩阵 (Nx × Ny)，中心化傅里叶变换结果，复数
%
% 【参考文献】
%   Tatiana Latychevskaia and Hans-Werner Fink
%   "Practical algorithms for simulation and reconstruction of digital in-line holograms",
%   Appl. Optics 54, 2424-2434 (2015)
%
% 【作者】Tatiana Latychevskaia, 2002 (MATLAB R2010b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = FT2Dc(in)

% 获取输入矩阵的行列数
[Nx, Ny] = size(in);

% 构造行索引列向量 ii (Nx×1) 和列索引行向量 jj (1×Ny)
% cast(..., 'like', in) 保证与输入相同的数据类型（如single/double）
ii = cast((1:Nx)', 'like', in);
jj = cast(1:Ny, 'like', in);

% 中心化移位因子 f1 = exp(jπ(m+n)) = (-1)^(m+n)
% 利用 MATLAB 广播机制，ii+jj 生成 Nx×Ny 矩阵
f1 = exp(1i * pi * (ii + jj));

% 步骤1: 输入乘以移位因子 (-1)^(m+n)，将零频移至频谱中心
% 步骤2: 执行标准二维FFT
FT = fft2(f1.*in);

% 步骤3: 输出再乘以移位因子 (-1)^(m+n)，完成空间域的中心化对齐
out = f1.*FT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
