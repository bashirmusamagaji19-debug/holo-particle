function [Phase3D Pupil]=MyMakingPhase3D(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,sensor_size)
% MyMakingPhase3D  为所有深度层批量生成频域传播相位和光瞳函数（3D体）。
%
% 【在CS链路中的位置】
%   在AS_reconstructCS中，于构造传播算子之前被调用，一次性生成全部层的
%   Phase3D和Pupil，供后续 MyFieldsPropagation / MyForwardPropagation /
%   MyAdjointPropagation 使用。
%   调用链：AS_reconstructCS → MyMakingPhase3D → (循环)MyMakingPhase
%
% 【输入】
%   Nx          - 频域网格行数（补零后全息图行数）
%   Ny          - 频域网格列数
%   Nz          - 深度层数
%   lambda      - 波长 [μm]
%   deltaX      - 空间域x方向采样间隔 [μm]
%   deltaY      - 空间域y方向采样间隔 [μm]
%   deltaZ      - 层间深度步长 [μm]
%   offsetZ     - 首层距探测器的偏移距离 [μm]（zMin对应的微米值）
%   sensor_size - 传感器物理尺寸 [μm]，用于计算各层NA
%
% 【输出】
%   Phase3D - Nx×Ny×Nz 复数3D数组，每层(:,:,i)为该深度层的频域传播相位
%   Pupil   - Nx×Ny×Nz 实数3D数组，每层(:,:,i)为该层的光瞳掩模
%
% 【物理说明】
%   NA随深度变化：NA = (sensor_size/2) / sqrt(z² + (sensor_size/2)²)
%   距离越远，NA越小，对应更大的衍射受限角

% 构建深度轴 Z = [offsetZ, offsetZ+deltaZ, ..., offsetZ+(Nz-1)*deltaZ]
Z=[0:Nz-1].*deltaZ+offsetZ;

% 传感器半尺寸，用于计算各层NA
s=sensor_size/2;

% 预分配3D输出数组
Phase3D=zeros(Nx,Ny,Nz);
Pupil=zeros(Nx,Ny,Nz);

% 逐层调用 MyMakingPhase 生成相位和光瞳
for i=1:Nz
    NA=s/sqrt(Z(i).^2+s.^2);  % 根据层深计算该层的数值孔径
    [Phase3D(:,:,i) Pupil(:,:,i)]=...
        MyMakingPhase(Nx,Ny,Z(i),lambda,deltaX,deltaY,NA);
end;
