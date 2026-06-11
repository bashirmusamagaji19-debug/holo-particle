function S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil)
% MyForwardPropagation  正向传播模型：从3D散射势eta预测探测器上的全息图S。
%
% 【在CS链路中的位置】
%   实现CS正问题 S = A(eta) 的核心物理计算。
%   被MyForwardOperatorPropagation包装后作为TwIST的正向算子A传入。
%   调用链：TwIST迭代 → A(f_twist) → MyForwardOperatorPropagation → 本函数
%
% 【物理原理】基于一阶Born近似的正向散射模型：
%   1. 散射场 = 散射势 × 照明场：Es = eta .* E
%   2. 各层散射场传播到探测器：FFT(Es) · Phase3D · Pupil
%   3. 所有层的贡献相加（相干叠加）：沿第3维求和
%   4. IFFT回空间域，取实部得到探测器强度 S
%
% 【输入】
%   eta     - Nx×Ny×Nz 复数3D散射势（待重建的未知量）
%   E       - Nx×Ny×Nz 复数3D照明场（由MyFieldsPropagation生成）
%   Nx,Ny   - 空间域网格尺寸
%   Nz      - 深度层数
%   phase3D - Nx×Ny×Nz 频域传播相位3D体
%   pupil   - Nx×Ny×Nz 光瞳函数3D体
%
% 【输出】
%   S       - Nx×Ny 实数2D矩阵，模拟的探测器全息图

cEs=zeros(Nx,Ny,Nz,'like',eta);
Es=eta.*E;  % 散射场 = 散射势 × 照明场（逐元素相乘，Born一阶近似）

% 各层散射场变换到频域
for i=1:Nz
    cEs(:,:,i)=fftshift(fft2(Es(:,:,i)));
end

% 频域中：各层传播到探测器（乘以Phase3D和Pupil），然后相干叠加（沿深度求和）
cEsp=sum(cEs.*phase3D.*pupil,3);

% IFFT回空间域，取实部作为探测器测量值
S=real((ifft2(ifftshift(cEsp))));
