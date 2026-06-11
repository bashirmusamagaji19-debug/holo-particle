function y=MyC2V(x)
% MyC2V  将复数向量展成实数向量，供TwIST等只接受实数输入的优化器使用。
%
% 【在CS链路中的位置】
%   TwIST只能处理实向量，而全息重建中的物理量（散射势 eta、观测全息图 g）都是复数。
%   MyC2V 把复数向量拆成 [实部; 虚部] 的实向量，使数据能在TwIST框架内传递。
%   调用链：AS_reconstructCS → 构造g_vec=MyC2V(g(:)) → TwIST → MyV2C恢复复数
%
% 【输入】
%   x - 复数列向量，长度为 N
%
% 【输出】
%   y - 实数列向量，长度为 2N，格式为 [real(x); imag(x)]
%
% 【数据流示意】
%   复数 z = a + bi  →  实向量 [a; b]
%   恢复时使用 MyV2C 做逆变换

y=[real(x);imag(x)];  % 上半段存放实部，下半段存放虚部
