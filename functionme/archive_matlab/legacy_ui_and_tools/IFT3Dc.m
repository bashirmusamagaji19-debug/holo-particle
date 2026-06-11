%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 三维中心化逆傅里叶变换 (3D Centered Inverse Fourier Transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 【用途】
%   对三维中心化频谱执行逆傅里叶变换，恢复空间域三维信号。
%   是 FT3Dc 的逆操作，用于三维全息重建中从频域恢复空域体数据。
%
% 【数学公式】
%   三维中心化逆变换：
%     g(m,n,k) = (-1)^(m+n+k) · IFFT3{ (-1)^(m+n+k) · G(u,v,w) }
%   移位因子指数取负：exp(-jπ(m+n+k))
%
% 【输入】
%   in  - 三维数组 (Nx × Ny × Nz)，中心化频谱，复数
%
% 【输出】
%   out - 三维数组 (Nx × Ny × Nz)，空间域复振幅
%
% 【参考文献】
%   Tatiana Latychevskaia and Hans-Werner Fink
%   "Practical algorithms for simulation and reconstruction of digital in-line holograms",
%   Appl. Optics 54, 2424-2434 (2015)
%
% 【作者】Tatiana Latychevskaia, 2002 (MATLAB R2010b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = IFT3Dc(in)

% 获取三维数组的各维尺寸
[Nx Ny Nz] = size(in);

% 预分配逆移位因子数组
f1 = zeros(Nx,Ny,Nz);

% 三重循环构造逆移位因子 f1(ii,jj,kk) = exp(-jπ(ii+jj+kk)) = (-1)^(ii+jj+kk)
% 注意：指数符号为负，与正变换 FT3Dc 相反
for ii = 1:Nx
    for jj = 1:Ny
        for kk = 1:Nz
        f1(ii, jj, kk) = exp(-i*pi*(ii + jj + kk));
        end
    end
end

% 步骤1: 输入频谱乘以逆移位因子
% 步骤2: 执行三维逆FFT (ifftn)
FT = ifftn(f1.*in);

% 步骤3: 输出再乘以逆移位因子，完成中心化恢复
out = f1.*FT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
