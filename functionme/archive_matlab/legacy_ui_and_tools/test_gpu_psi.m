%TEST_GPU_PSI GPU上MyTVpsi函数的测试脚本
%
%   【用途】
%   测试MyTVpsi函数在GPU上的运行情况。
%   MyTVpsi是全变分（TV）正则化中的梯度算子，用于压缩感知迭代重建。
%   本脚本验证该函数能否正确接收GPU数组输入并返回正确尺寸的输出。
%
%   【在工程中的位置】
%   独立测试脚本，不在主链路中。用于开发调试阶段验证GPU兼容性。
%
%   【注意】
%   首行cd路径包含中文字符和原始开发环境路径，需根据实际环境修改。

% 切换到函数目录（路径需按实际环境修改）
cd(''D:/周柏臻/UI界面To苏老师/function'');
try
    % 定义小尺寸测试参数
    nx = 8; ny = 8; nz = 3;
    % 在GPU上生成随机单精度测试向量（长度=2*nx*ny*nz，对应实部+虚部展开）
    x = gpuArray.rand(2*nx*ny*nz,1,'single');
    % 调用MyTVpsi：输入x，正则化参数0.01，权重0.05，维度参数1, nx, ny*nz*2, 1
    y = MyTVpsi(x, single(0.01), 0.05, 1, nx, ny*nz*2, 1);
    % 打印输出尺寸，验证正确性
    disp(size(y));
catch ME
    % 捕获并显示错误详情
    disp(getReport(ME,'extended','hyperlinks','off'));
end
