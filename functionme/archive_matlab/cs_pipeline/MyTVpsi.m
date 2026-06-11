function y=MyTVpsi(x,th,tau,iter,Nvx,Nvy,Nvz)
% MyTVpsi - TV正则项的近端收缩接口 Psi(x, th)
%
% 【用途】
%   计算 TV 正则项的 proximal 操作（近端映射/软收缩）：
%     Ψ(x, θ) = prox_{θ·TV}(x) = argmin_z { (1/2)||z-x||² + θ·TV(z) }
%   这是 TwIST 迭代中用于正则项更新的核心操作。
%
% 【在TV正则化链路中的位置】
%   被 TwIST 优化器直接调用，作为 Psi 接口。
%   调用链：MyTVpsi → MyProjectionTV → MyTV3D_conv / MyDiv3D
%
% 【数学含义】
%   近端算子的对偶表示（Chambolle 分解）：
%     prox_{θ·TV}(x) = x - θ·div(p*)
%   其中 p* 是对偶最优解，通过 MyProjectionTV 迭代求解。
%   等价地：Y = X - MyProjectionTV(X, τ, θ/2, iter)
%   注意：传入的 λ = θ/2，这是因为目标函数的系数关系。
%
%   直观理解：proximal 操作在保持与 x 接近的同时，
%   使输出更"分段常数"（TV正则的先验效果）。
%
% 【输入】
%   x    - 列向量 (Nvx*Nvy*Nvz, 1)，当前迭代解（一维存储）
%   th   - 标量，正则化参数 θ（阈值），控制收缩强度
%          θ 越大，收缩越强，结果越平滑/稀疏
%   tau  - 标量，对偶迭代的步长参数
%   iter - 正整数，对偶迭代次数
%   Nvx  - x 方向维度大小
%   Nvy  - y 方向维度大小
%   Nvz  - z 方向维度大小
%
% 【输出】
%   y - 列向量 (Nvx*Nvy*Nvz, 1)，近端收缩后的结果（一维存储）

% 将一维向量重塑为三维数组
X=reshape(x,Nvx,Nvy,Nvz);

% 计算 proximal 操作：Y = X - ProjectionTV(X, τ, λ, iter)
% MyProjectionTV 返回 λ·div(p*)，即被收缩掉的分量
% 因此 Y = X - λ·div(p*) = prox_{θ·TV}(X)
% 传入 lam = th*0.5 是因为 proximal 目标函数中系数为 θ/2
Y=X-MyProjectionTV(X,tau,th*0.5,iter);

% 将三维结果重塑回一维列向量，返回给 TwIST 优化器
y=reshape(Y,Nvx*Nvy*Nvz,1);
