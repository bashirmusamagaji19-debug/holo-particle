function p=MyProjectionTV(g,tau,lam,iter)
% MyProjectionTV - TV对偶投影 / Proximal算子内核
%
% 【用途】
%   通过对偶方法求解 TV 的 proximal 操作（近端映射）：
%     prox_{λ·TV}(g) = g - λ · div(p*)
%   其中 p* 是以下对偶问题的解：
%     p* = argmin_p  ||g - λ·div(p)||²  s.t. ||p||_∞ ≤ 1
%
% 【在TV正则化链路中的位置】
%   MyTVpsi → MyProjectionTV  （被 MyTVpsi 调用，计算近端收缩结果）
%   内部调用 MyTV3D_conv 和 MyDiv3D 进行梯度/散度计算。
%
% 【数学含义】
%   采用 Chambolle 对偶迭代法求解 TV proximal 问题。
%   迭代格式（梯度上升法 + 投影）：
%     a = ∇(div(p_n) - g/λ)
%     p_{n+1} = (p_n + τ·a) / (1 + τ·||a||_2)
%   其中分母 1+τ·||a||_2 实现向单位球 {p: ||p||_∞≤1} 的投影。
%
%   算法参考：Chambolle, "An Algorithm for Total Variation
%   Minimization and Applications", JMIV 2004.
%
% 【输入】
%   g    - 三维数组 (nx, ny, nz)，输入图像（待收缩的信号）
%   tau  - 标量，对偶迭代的步长参数（需 τ < 1/8 保证收敛，3D情况）
%   lam  - 标量，正则化参数 λ（在 MyTVpsi 中传入 th*0.5）
%   iter - 正整数，对偶迭代次数
%
% 【输出】
%   p - 三维数组 (nx, ny, nz)，proximal 结果 = g - λ·div(p*)
%       即 g 减去 λ 倍的对偶变量散度

[nx,ny,nz]=size(g);

% 初始化对偶变量 p_n = 0（三维向量场，3个分量）
pn=zeros(nx,ny,nz,3,'like',g);
% div(p_n) 的缓存，初始为零
div_pn=zeros(nx,ny,nz,'like',g);
% b 用于存储梯度幅值，辅助投影计算
b=pn;

% ---- Chambolle 对偶迭代 ----
% 迭代求解约束优化：max_p <g, div(p)> - (λ/2)||div(p)||²  s.t. ||p||_∞ ≤ 1
for i=1:iter

    % 计算目标函数的梯度方向：
    % a = ∇(div(p_n) - g/λ) = ∇(div_pn - g/lam)
    % 这是上升方向，用于更新对偶变量
    a=MyTV3D_conv(div_pn-g./lam);

    % 计算梯度幅值 ||a||_2 = √(ax² + ay² + az²)
    % 将标量幅值扩展到3个分量，用于逐分量投影
    b(:,:,:,1)=sqrt(a(:,:,:,1).^2+a(:,:,:,2).^2+a(:,:,:,3).^2);
    b(:,:,:,2)=b(:,:,:,1);
    b(:,:,:,3)=b(:,:,:,1);

    % 投影梯度上升更新：
    % p_{n+1} = (p_n + τ·a) / (1 + τ·||a||_2)
    % 分母确保 p 被投影到单位球 {p: ||p||_2 ≤ 1} 内
    % 这是逐点投影，每个体素独立
    pn=(pn+tau.*a)./(1.0+tau.*b);

    % 更新散度 div(p_n)，供下一次迭代使用
    div_pn=MyDiv3D(pn);
end;

% 输出 proximal 结果：g - λ·div(p*)
% 其中 λ·div(p*) 是被收缩掉的部分（TV惩罚对应的信号分量）
p=lam.*MyDiv3D(pn);
