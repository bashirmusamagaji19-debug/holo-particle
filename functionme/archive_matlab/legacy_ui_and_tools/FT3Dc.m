%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 三维中心化傅里叶变换 (3D Centered Fourier Transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 【用途】
%   对三维输入数组执行中心化傅里叶变换，使零频分量位于三维频谱中心。
%   用于三维全息重建中体数据的频域分析和滤波。
%
% 【数学公式】
%   三维中心化移位性质：
%     G(u,v,w) = (-1)^(m+n+k) · FFT3{ (-1)^(m+n+k) · g(m,n,k) }
%   其中 (-1)^(m+n+k) = exp(jπ(m+n+k))，原理与二维情形相同。
%
% 【输入】
%   in  - 三维数组 (Nx × Ny × Nz)，可为实数或复数
%
% 【输出】
%   out - 三维数组 (Nx × Ny × Nz)，中心化三维傅里叶变换结果，复数
%
% 【参考文献】
%   Tatiana Latychevskaia and Hans-Werner Fink
%   "Practical algorithms for simulation and reconstruction of digital in-line holograms",
%   Appl. Optics 54, 2424-2434 (2015)
%
% 【作者】Tatiana Latychevskaia, 2002 (MATLAB R2010b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = FT3Dc(in)

% 获取三维数组的各维尺寸
[Nx, Ny, Nz] = size(in);

% 预分配移位因子数组
f1 = zeros(Nx,Ny,Nz);

% 三重循环构造移位因子 f1(ix,iy,iz) = exp(jπ(ix+iy+iz)) = (-1)^(ix+iy+iz)
% 注：此处索引从1开始，(-1)^(ix+iy+iz) 仍满足移位性质
for ix = 1:Nx
    for iy = 1:Ny
          for iz = 1:Nz
              f1(ix,iy,iz) = exp(1i*pi*(ix + iy + iz));
          end
    end
end

% 步骤1: 输入乘以移位因子
% 步骤2: 执行三维FFT (fftn)
FT = fftn(f1.*in);

% 步骤3: 输出再乘以移位因子，完成中心化
out = f1.*FT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
