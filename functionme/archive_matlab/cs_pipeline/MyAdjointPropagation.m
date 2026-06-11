function eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil)
% MyAdjointPropagation  伴随传播模型：正向算子A的共轭转置AT的实现。
%
% 【在CS链路中的位置】
%   实现CS逆问题中梯度计算所需的伴随算子 AT。
%   被MyAdjointOperatorPropagation包装后作为TwIST的AT参数传入。
%   TwIST每次迭代计算梯度：grad = AT(resid) = AT(y - A(x))
%   调用链：TwIST迭代 → AT(gmeas) → MyAdjointOperatorPropagation → 本函数
%
% 【物理/数学原理】
%   正向模型：S = Re{ IFFT{ Σ_z FFT{eta_z·E_z} · Phase3D_z · Pupil_z } }
%   伴随算子是正向算子的共轭转置，数学上等价于将正向流程"逆序+共轭"：
%   1. 将探测器全息图S变换到频域
%   2. 频域中广播到所有层，乘以共轭相位和共轭光瞳
%   3. 各层IFFT回空间域
%   4. 乘以照明场的共轭：eta = conj(E) .* result
%
% 【输入】
%   S       - Nx×Ny 实数2D矩阵，探测器测量（或残差）
%   E       - Nx×Ny×Nz 复数3D照明场
%   Nx,Ny   - 空间域网格尺寸
%   Nz      - 深度层数
%   phase3D - Nx×Ny×Nz 频域传播相位3D体
%   pupil   - Nx×Ny×Nz 光瞳函数3D体
%
% 【输出】
%   eta     - Nx×Ny×Nz 复数3D散射势的伴随重建（梯度方向）

% 步骤1：将探测器信号S变换到频域
% 先对S取实部、取共轭、做IFFT（等效于正向IFFT的共轭转置），再fftshift
cEsp=fftshift(conj(ifft2(conj(real(S)))));

% 步骤2：频域中，各层乘以共轭传播相位和共轭光瞳（正向乘Phase*Pupil的共轭转置）
cEs=conj(phase3D).*conj(pupil).*repmat(cEsp,[1 1 Nz]);  % 广播到Nz层

% 步骤3：各层从频域变换回空间域（正向FFT的共轭转置）
eta=zeros(Nx,Ny,Nz,'like',cEs);
for i=1:Nz
    eta(:,:,i)=conj(fft2(conj(ifftshift(cEs(:,:,i)))));
end

% 步骤4：乘以照明场的共轭（正向中 eta.*E 的共轭转置为 conj(E).*result）
eta=conj(E).*eta;
