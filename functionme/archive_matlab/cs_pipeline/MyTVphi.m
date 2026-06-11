function y=MyTVphi(x,Nvx,Nvy,Nvz)
% MyTVphi - TV正则项接口函数 Phi(x)
%
% 【用途】
%   计算 TV 正则项的值：Φ(x) = TV(x) = Σ ||∇x||_2
%   该值用于 TwIST 优化器的停止准则判断
%   （当 Φ(x) 变化量小于阈值时，认为算法收敛）。
%
% 【在TV正则化链路中的位置】
%   被 TwIST 优化器直接调用，作为 Phi 接口。
%   调用链：MyTVphi → MyTVnorm → MyTV3D_conv
%
% 【数学含义】
%   Φ(x) = Σ_{i,j,k} √( (∂x/∂i)² + (∂x/∂j)² + (∂x/∂k)² )
%   即三维各向同性全变分范数。
%
% 【输入】
%   x    - 列向量 (Nvx*Nvy*Nvz, 1)，TwIST 中变量以一维向量形式存储
%   Nvx  - x 方向维度大小
%   Nvy  - y 方向维度大小
%   Nvz  - z 方向维度大小
%
% 【输出】
%   y - 标量，TV 正则项值 Φ(x)

% 将一维向量重塑为三维数组，以适配三维差分运算
X=reshape(x,Nvx,Nvy,Nvz);

% 计算三维TV范数，y 为标量总和，dif 为每个体素的梯度幅值
[y,dif]=MyTVnorm(X);
