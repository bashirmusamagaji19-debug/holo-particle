function y=MyDiv3D(TV)
% MyDiv3D - 三维散度算子
%
% 【用途】
%   计算三维向量场 TV 的散度，即 div(TV)。
%   这是 MyTV3D_conv（前向差分/梯度）的伴随算子（adjoint operator）。
%
% 【在TV正则化链路中的位置】
%   MyTVpsi → MyProjectionTV → MyDiv3D  （对偶迭代中计算散度）
%
%   与 MyTV3D_conv 满足伴随关系：<∇x, p> = <x, div p>，
%   其中 ∇ = MyTV3D_conv，div = MyDiv3D。
%
% 【数学含义】
%   散度是前向差分的伴随（后向差分）：
%     (div p)_{i,j,k} = p^x_{i,j,k} - p^x_{i-1,j,k}
%                      + p^y_{i,j,k} - p^y_{i,j-1,k}
%                      + p^z_{i,j,k} - p^z_{i,j,k-1}
%
%   边界条件（与 MyTV3D_conv 的 Neumann 边界配对）：
%     - 第1个像素：无 i-1 邻居，只取正向项 p^x_{1,j,k}
%     - 最后1个像素：取负向项 -p^x_{N-1,j,k}（循环偏移后的值）
%     - y/z 方向同理
%
% 【输入】
%   TV - 四维数组 (nx, ny, nz, 3)，三维向量场
%        TV(:,:,:,1) = x 分量
%        TV(:,:,:,2) = y 分量
%        TV(:,:,:,3) = z 分量
%
% 【输出】
%   y - 三维数组 (nx, ny, nz)，散度场 div(TV)

n=size(TV);

% ---- x 方向散度：∂p^x/∂i ≈ p^x_{i,j,k} - p^x_{i-1,j,k} ----
% circshift(TV(:,:,:,1),[1,0,0]) 等价于 p^x_{i-1,j,k}（后向偏移）
x_shift=circshift(TV(:,:,:,1),[1 0 0]);
yx=TV(:,:,:,1)-x_shift;
% 边界 i=1：无前向邻居，散度 = p^x_{1,j,k}
yx(1,:,:)=TV(1,:,:,1);
% 边界 i=nx：取 -p^x_{nx-1,j,k}（循环偏移值修正）
yx(n(1),:,:)=-x_shift(n(1),:,:);

% ---- y 方向散度：∂p^y/∂j ≈ p^y_{i,j,k} - p^y_{i,j-1,k} ----
y_shift=circshift(TV(:,:,:,2),[0 1 0]);
yy=TV(:,:,:,2)-y_shift;
% 边界 j=1
yy(:,1,:)=TV(:,1,:,2);
% 边界 j=ny
yy(:,n(2),:)=-y_shift(:,n(2),:);

% ---- z 方向散度：∂p^z/∂k ≈ p^z_{i,j,k} - p^z_{i,j,k-1} ----
z_shift=circshift(TV(:,:,:,3),[0 0 1]);
yz=TV(:,:,:,3)-z_shift;
% 边界 k=1
yz(:,:,1)=TV(:,:,1,3);
% 边界 k=nz
yz(:,:,n(3))=-z_shift(:,:,n(3));

% ---- 总散度 = 三个方向散度之和 ----
y=yx+yy+yz;
