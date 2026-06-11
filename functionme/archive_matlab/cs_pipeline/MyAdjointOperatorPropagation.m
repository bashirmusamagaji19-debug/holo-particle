function eta=MyAdjointOperatorPropagation(S,E,Nx,Ny,Nz,phase3D,pupil)
% MyAdjointOperatorPropagation  伴随传播算子的向量化接口，适配TwIST函数句柄调用。
%
% 【在CS链路中的位置】
%   作为TwIST的伴随算子AT的函数句柄实现。TwIST内部计算梯度时调用
%   AT(gmeas)，执行此函数，负责将TwIST的实数向量输入转换为2D复数矩阵，
%   调用物理伴随传播，再将3D输出展平为实数向量返回。
%   调用链：TwIST → AT(gmeas) → 本函数 → MyAdjointPropagation
%
% 【数据转换流程】
%   TwIST实向量 gmeas [2·Nx·Ny × 1]
%     → MyV2C: 实向量还原为复数向量 [Nx·Ny × 1]
%     → reshape: 重塑为2D复数矩阵 [Nx × Ny]
%     → MyAdjointPropagation: 物理伴随传播，输出3D散射势 [Nx × Ny × Nz]
%     → 展平 → MyC2V: 复数展为实向量 [2·Nx·Ny·Nz × 1] 返回给TwIST
%
% 【输入】
%   S       - 实数列向量 [2·Nx·Ny × 1]，TwIST的残差/测量数据
%             格式为 [real; imag]（由MyC2V编码）
%   E       - Nx×Ny×Nz 复数3D照明场
%   Nx,Ny   - 空间域网格尺寸
%   Nz      - 深度层数
%   phase3D - Nx×Ny×Nz 频域传播相位3D体
%   pupil   - Nx×Ny×Nz 光瞳函数3D体
%
% 【输出】
%   eta     - 实数列向量 [2·Nx·Ny·Nz × 1]，伴随传播结果的[real;imag]编码

% 实向量 → 复数向量 → 2D复数矩阵
S=reshape(MyV2C(S),Nx,Ny);

% 调用物理伴随传播
eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil);

% 3D散射势 → 列向量 → 实向量编码
eta=MyC2V(eta(:));
