function S=MyForwardOperatorPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil)
% MyForwardOperatorPropagation  正向传播算子的向量化接口，适配TwIST函数句柄调用。
%
% 【在CS链路中的位置】
%   作为TwIST的正向算子A的函数句柄实现。TwIST内部调用 A(f_twist) 时
%   执行此函数，负责将TwIST的实数向量输入转换为3D复数体，调用物理正向
%   传播，再将输出展平为实数向量返回。
%   调用链：TwIST → A(f_twist) → 本函数 → MyForwardPropagation
%
% 【数据转换流程】
%   TwIST实向量 f_twist [2·Nx·Ny·Nz × 1]
%     → MyV2C: 实向量还原为复数向量 [Nx·Ny·Nz × 1]
%     → reshape: 重塑为3D复数体 [Nx × Ny × Nz]
%     → MyForwardPropagation: 物理正向传播，输出2D全息图 [Nx × Ny]
%     → MyC2V: 复数展为实向量 [2·Nx·Ny × 1] 返回给TwIST
%
% 【输入】
%   eta     - 实数列向量 [2·Nx·Ny·Nz × 1]，TwIST的当前估计
%             格式为 [real; imag]（由MyC2V编码）
%   E       - Nx×Ny×Nz 复数3D照明场
%   Nx,Ny   - 空间域网格尺寸
%   Nz      - 深度层数
%   phase3D - Nx×Ny×Nz 频域传播相位3D体
%   pupil   - Nx×Ny×Nz 光瞳函数3D体
%
% 【输出】
%   S       - 实数列向量 [2·Nx·Ny × 1]，正向传播结果的[real;imag]编码

% 实向量 → 复数向量 → 3D复数体
eta=reshape(MyV2C(eta),Nx,Ny,Nz);

% 调用物理正向传播
S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil);

% 2D全息图 → 列向量 → 实向量编码
S=MyC2V(S(:));
