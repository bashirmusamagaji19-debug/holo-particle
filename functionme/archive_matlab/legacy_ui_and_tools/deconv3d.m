%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 三维反卷积 (3D Deconvolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 【用途】
%   从全息强度重建结果中去除点扩展函数(PSF)的模糊效应，
%   提高三维重建的轴向和横向分辨率。采用维纳滤波(Wiener filter)
%   形式的频域反卷积。
%
% 【数学公式】
%   全息重建的强度像可近似为：
%     I_3D(r) ≈ |U_obj(r)|² ⊗ |PSF(r)|²
%   频域中：
%     Î(fx,fy,fz) = Ô(fx,fy,fz) · CTF(fx,fy,fz)
%   其中 CTF = FT{|PSF|²} 为对比度传递函数。
%
%   维纳反卷积滤波器：
%     CTFi(fx,fy,fz) = CTF*(fx,fy,fz) / (|CTF(fx,fy,fz)|² + β)
%   β 为正则化参数，防止 CTF 零点处除零并抑制噪声放大。
%
%   反卷积结果：
%     obj_deconv = IFT{ FT{I_3D} · CTFi }
%
%   对应文献 Eq.(20)：Opt. Express 18(21):22527-22544, 2010
%
% 【适用条件】
%   当各散射点之间的干涉效应可忽略时成立，例如粒子场全息中
%   大量相同粒子分散在体积内的情形。
%
% 【输入】
%   obj_3d_complex - 三维复振幅重建结果 (Nx × Ny × Nz)，复数
%   PSF            - 三维点扩展函数 (Nx × Ny × Nz)，复数
%   beta           - 正则化参数（小正数），避免除零和抑制噪声
%
% 【输出】
%   obj_3d_deconv  - 三维反卷积结果 (Nx × Ny × Nz)
%
% 【参考文献】
%   Section 4.1 in Opt. Express 18(21):22527-22544, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obj_3d_deconv]= deconv3d(obj_3d_complex, PSF, beta)

    % 定义三维中心化FFT和逆FFT
    FTn = @(x) ifftshift(fftn(fftshift(x)));
    iFTn = @(x) ifftshift(ifftn(fftshift(x)));
    
    % 取复振幅的强度（模的平方）作为重建强度像
    % I_3D(r) = |U(r)|²
    obj_3d = abs(obj_3d_complex).^2;

    % 计算 CTF = FT{|PSF|²}，即PSF强度的傅里叶变换
    CTF = FTn(abs(PSF).^2);
    clear PSF

    % 构造维纳反卷积滤波器
    % CTFi = CTF* / (|CTF|² + β)
    % β 是小正数：避免 CTF 零点处除零，同时抑制噪声放大
    % β 越大，正则化越强，结果越平滑但细节越少
    CTFi = conj(CTF)./(abs(CTF).^2 + beta);
    clear CTF

    % 执行反卷积：频域相乘后逆变换
    % obj_3d_deconv = IFT{ FT{I_3D} · CTFi }  — 文献 Eq.(20)
    obj_3d_deconv = iFTn(FTn(obj_3d).*CTFi);
    clear CTFi
end
