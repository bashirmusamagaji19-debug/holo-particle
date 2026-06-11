%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rayleigh-Sommerfeld 角谱传播法 (Angular Spectrum Method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
【用途】
  基于角谱理论的自由空间波前传播，将输入光场 U 沿 z 方向传播距离 z。
  该方法为严格标量衍射理论（Rayleigh-Sommerfeld），不做傍轴近似，
  适用于任意传播距离（近场和远场），但要求采样满足奈奎斯特条件。

【数学公式】
  角谱传播的传递函数：
    P(fx, fy) = exp(-j · 2π · z / λ · sqrt(1 - λ²(fx² + fy²)))
  其中：
    fx, fy  — 空间频率
    λ       — 波长
    z       — 传播距离
    α = λ·fx, β = λ·fy — 归一化方向余弦

  当 α² + β² > 1 时，对应倏逝波，传递函数指数衰减（虚数开根号产生衰减）。

  传播过程：U_z = IFT{ FT{U} · P }

【输入】
  U     - 输入复振幅光场 (N × N)，中心化排列
  N     - 采样点数（假设 x 和 y 方向相同）
  lambda - 波长 (m)
  area  - 计算窗口物理尺寸 (m)，即 N × dx
  z     - 传播距离 (m)

【输出】
  U_z   - 传播后的复振幅光场 (N × N)

【参考文献】
  Tatiana Latychevskaia and Hans-Werner Fink
  "Practical algorithms for simulation and reconstruction of digital in-line holograms",
  Appl. Optics 54, 2424-2434 (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

function [U_z] = RS(U, N, lambda, area, z)

    % 定义中心化FFT和逆FFT的快捷函数
    % fftshift 将零频移至中心，ifftshift 为其逆操作
    % 组合使用：ifftshift(fft2(fftshift(x))) 实现中心化FFT
    FT2 = @(x) ifftshift(fft2(fftshift(x)));
    iFT2 = @(x) ifftshift(ifft2(fftshift(x)));

    % 构造空间频率坐标网格（以零频为中心）
    % (1:N)-N/2-1 生成 [-N/2, N/2-1] 的整数坐标
    [x, y] = meshgrid((1:N)-N/2-1);

    % 计算归一化方向余弦
    % alpha = λ·fx = λ·(像素索引)/area，beta 同理
    alpha = lambda.*x./area;
    beta = lambda.*y./area;

    % 计算 r = α² + β²，用于判断传播波和倏逝波
    r = alpha.^2 + beta.^2;

    % 角谱传递函数
    % 当 r ≤ 1 时，sqrt(1-r) 为实数，对应传播波（相位调制）
    % 当 r > 1 时，sqrt(1-r) 为虚数，对应倏逝波（指数衰减）
    p = exp(-2*pi*1i*z.*sqrt(1-r)./lambda);

    % 倏逝波处理：可取消下行注释将倏逝波分量置零
    % p(isnan(p)) = 0;

    % 传播计算：频域相乘，再逆变换回空域
    U_z = iFT2(FT2(U).*p);
end
