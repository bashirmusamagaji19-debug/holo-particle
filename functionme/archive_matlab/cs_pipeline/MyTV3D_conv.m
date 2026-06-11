function TV=MyTV3D_conv(x)
% MyTV3D_conv - 三维前向差分算子（离散梯度）
%
% 【用途】
%   计算三维图像 x 沿 x/y/z 三个方向的前向差分，即离散梯度 ∇x。
%   结果 TV(:,:,:,d) 存储第 d 个方向（1=x, 2=y, 3=z）的差分。
%
% 【在TV正则化链路中的位置】
%   MyTVphi → MyTVnorm → MyTV3D_conv   （计算梯度，用于求TV范数）
%   MyTVpsi → MyProjectionTV → MyTV3D_conv （对偶迭代中计算梯度）
%
%   与 MyDiv3D 互为伴随算子（adjoint），即 <∇x, p> = <x, div p>，
%   这是 TV 对偶方法的核心性质。
%
% 【数学含义】
%   前向差分定义（Neumann 边界条件）：
%     (∇x)_{i,j,k,1} = x_{i+1,j,k} - x_{i,j,k}   （x方向）
%     (∇x)_{i,j,k,2} = x_{i,j+1,k} - x_{i,j,k}   （y方向）
%     (∇x)_{i,j,k,3} = x_{i,j,k+1} - x_{i,j,k}   （z方向）
%   边界处（最后一个像素）差分置零，对应 Neumann 边界条件 ∂x/∂n = 0。
%
% 【输入】
%   x - 三维数组 (nx, ny, nz)，待求梯度的图像/信号
%
% 【输出】
%   TV - 四维数组 (nx, ny, nz, 3)，第4维为三个方向的差分
%        TV(:,:,:,1) = x方向前向差分
%        TV(:,:,:,2) = y方向前向差分
%        TV(:,:,:,3) = z方向前向差分

[nx,ny,nz]=size(x);

% 初始化四维输出数组，与输入同类型（支持 GPU/single/double）
TV=zeros(nx,ny,nz,3,'like',x);

% ---- x 方向前向差分 ----
% circshift(x,[-1,0,0]) 将 x 沿第1维向上平移1，等价于 x_{i+1,j,k}
% 差分 = x_{i+1,j,k} - x_{i,j,k}
TV(:,:,:,1)=circshift(x,[-1 0 0])-x;
% 边界处理：最后一行无 i+1 邻居，差分置零（Neumann边界）
TV(nx,:,:,1)=0.0;

% ---- y 方向前向差分 ----
% circshift(x,[0,-1,0]) 将 x 沿第2维向左平移1，等价于 x_{i,j+1,k}
TV(:,:,:,2)=circshift(x,[0 -1 0])-x;
% 边界处理：最后一列无 j+1 邻居，差分置零
TV(:,ny,:,2)=0.0;

% ---- z 方向前向差分 ----
% circshift(x,[0,0,-1]) 将 x 沿第3维平移1，等价于 x_{i,j,k+1}
TV(:,:,:,3)=circshift(x,[0 0 -1])-x;
% 边界处理：最后一层无 k+1 邻居，差分置零
TV(:,:,nz,3)=0.0;
% 乘以 1.0 保持数据类型一致（GPU数组类型兼容）
TV(:,:,:,3)=TV(:,:,:,3).*(1.0);
