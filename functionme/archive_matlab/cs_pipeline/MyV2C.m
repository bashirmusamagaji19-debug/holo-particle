function y=MyV2C(x)
% MyV2C  将实数向量还原为复数向量，是 MyC2V 的逆操作。
%
% 【在CS链路中的位置】
%   TwIST求解完成后输出实数向量 f_vec，需通过 MyV2C 将 [real;imag] 形式
%   还原为复数向量，再 reshape 成3D复数体并取幅值。
%   调用链：TwIST → f_c = MyV2C(f_vec) → reshape → abs → 3D幅值体
%
% 【输入】
%   x - 实数向量，长度必须为偶数 2N
%       前 N 个元素为复数实部，后 N 个元素为复数虚部
%
% 【输出】
%   y - 复数列向量，长度为 N，y = complex(x(1:N), x(N+1:2N))
%
% 【数据流示意】
%   实向量 [a; b]  →  复数 a + bi

x = x(:);                    % 确保为列向量
n = numel(x);
half = n / 2;
if rem(n, 2) ~= 0            % 长度必须为偶数，否则实虚部无法对半拆分
    error('MyV2C:InvalidLength', 'Input length must be even.');
end
y = complex(x(1:half), x(half+1:n));  % 前半段为实部，后半段为虚部
