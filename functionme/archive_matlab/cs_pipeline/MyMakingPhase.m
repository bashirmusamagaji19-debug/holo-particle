function [Phase Pupil]=MyMakingPhase(Nx,Ny,z,lambda,deltaX,deltaY,NA)
% MyMakingPhase  计算单层角谱传播相位传递函数与光瞳函数。
%
% 【在CS链路中的位置】
%   被MyMakingPhase3D逐层调用，为每一深度层z生成对应的频域相位因子和光瞳掩模。
%   这些相位和光瞳构成了正向/伴随传播算子的核心核函数。
%   调用链：AS_reconstructCS → MyMakingPhase3D → MyMakingPhase（每层一次）
%
% 【物理原理】
%   基于角谱法(Angular Spectrum Method)，在频域中：
%   - 传播相位：H(fx,fy;z) = exp(j·2π·z·sqrt(1/λ² - fx² - fy²))
%   - 光瞳函数：模拟有限NA的圆形低通滤波，限制传播角度范围
%
% 【输入】
%   Nx      - 频域网格行数（等于补零后的全息图行数）
%   Ny      - 频域网格列数
%   z       - 当前层传播距离 [μm]
%   lambda  - 波长 [μm]
%   deltaX  - 空间域x方向采样间隔 [μm]（补零后的有效像素尺寸）
%   deltaY  - 空间域y方向采样间隔 [μm]
%   NA      - 当前层对应的数值孔径
%
% 【输出】
%   Phase   - Nx×Ny 复数矩阵，频域传播相位因子 exp(j·2π·z·kz)
%   Pupil   - Nx×Ny 实数矩阵（目前为全1），可作为NA截断掩模使用

k=1/lambda;   % 波数 k = 1/λ（注意这里用了 1/λ 而非 2π/λ，后续相位公式中已包含2π）

% 构建频域坐标网格 fx, fy
X=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));  % fx 频率轴 [cycles/μm]
Y=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));   % fy 频率轴 [cycles/μm]

kx=repmat(X,1,Ny);   % kx(i,j) = fx(i)，频域x方向坐标网格
ky=repmat(Y,Nx,1);   % ky(i,j) = fy(j)，频域y方向坐标网格
kp=sqrt(kx.^2+ky.^2); % 径向频率 kp = sqrt(fx²+fy²)

% 计算传播项 kz = sqrt(k² - kp²)，其中 k=1/λ
% 当 kp > k 时对应倏逝波，令 term=0 予以滤除
term=k.^2-kp.^2;
term(term<0)=0;       % 倏逝波截止：负值置零，等效于低通滤波

% 频域传播相位因子 H = exp(j·2π·z·kz)
Phase=exp(j*2*pi*z*sqrt(term));

% 光瞳函数：当前实现为全1（无截断）
% 实际使用中可根据NA限制将 Pupil 设为圆形掩模
Pupil=ones(Nx,Ny);
