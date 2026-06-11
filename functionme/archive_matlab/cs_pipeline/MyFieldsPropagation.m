function E=MyFieldsPropagation(E0,Nx,Ny,Nz,phase3D,pupil)
% MyFieldsPropagation  计算照明光场E在各深度层的复振幅分布。
%
% 【在CS链路中的位置】
%   在AS_reconstructCS中，于构造传播算子之前被调用。
%   将均匀平面波E0(全1矩阵)通过角谱法传播到各深度层，得到照明场E。
%   E在正向/伴随算子中用于建模散射势与照明场的相互作用：S = A(eta) 中
%   正向传播先做 eta.*E 得到散射场，再传播回探测器。
%   调用链：AS_reconstructCS → MyFieldsPropagation → (循环)频域传播
%
% 【物理原理】
%   照明场E(z) = IFFT{ FFT{E0} · conj(Phase3D(:,:,z)) · Pupil(:,:,z) }
%   使用共轭相位 conj(Phase) 表示从探测器向深度的正向传播（向下传播）
%
% 【输入】
%   E0      - Nx×Ny 输入照明场（通常为全1平面波）
%   Nx,Ny   - 空间域网格尺寸
%   Nz      - 深度层数
%   phase3D - Nx×Ny×Nz 频域传播相位3D体（由MyMakingPhase3D生成）
%   pupil   - Nx×Ny×Nz 光瞳函数3D体
%
% 【输出】
%   E       - Nx×Ny×Nz 复数3D数组，照明场在各深度层的复振幅

E=zeros(Nx,Ny,Nz,'like',phase3D);
cE0=fftshift(fft2(E0));  % 对入射场做2D FFT并移中，得到频域表示

% 逐层计算照明场：频域乘以共轭传播相位和光瞳，再IFFT回空间域
for i=1:Nz
    % 共轭相位 conj(phase3D) 表示从源向深度的正向传播
    cE=cE0.*conj(phase3D(:,:,i)).*pupil(:,:,i);

    E(:,:,i)=ifft2(ifftshift(cE));  % 频域移回后再IFFT
end
