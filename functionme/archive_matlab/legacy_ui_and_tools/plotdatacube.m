function spatialcube = plotdatacube(data);
%PLOTDATACUBE 切片拼图函数
%
%   【用途】
%   将3D体数据的各层切片按网格排列拼接为一张2D大图，
%   便于一次性概览所有深度平面的重建结果。
%
%   【在工程中的位置】
%   辅助显示工具，供调试和快速查看3D体数据各层切片使用。
%   主链路不直接调用，但可在GUI的调试模式下使用。
%
%   【输入参数】
%       data : 3D体数据矩阵（ny × nx × nz），每层为一个切片
%
%   【输出参数】
%       spatialcube : 拼接后的2D大图（rows*ny × cols*nx），各切片按网格排列

    % 在每层切片四周填充3像素边界，填充值为数据最大值（使切片边界在显示中清晰可见）
    data = padarray(data,[3 3],max(data(:)),'both');

    n1 = size(data,1);  % 填充后切片行数
    n2 = size(data,2);  % 填充后切片列数

    totalfigs = size(data,3);  % 总切片数

    cols = 5;                           % 每行排列5个切片
    rows = ceil(totalfigs/cols);        % 计算所需行数

    % 预分配拼接后的大图矩阵
    spatialcube = zeros(rows*n1,cols*n2);

    figscount = 1;
    for r = 1:rows
        for c = 1:cols
            if figscount<=totalfigs
                % 将第figscount层切片放入大图对应位置
                spatialcube((r-1)*n1+1:(r-1)*n1+n1,(c-1)*n2+1:(c-1)*n2+n2) = squeeze(data(1:n1,1:n2,figscount));
            else
                % 切片数不足时用零填充空白位置
                spatialcube((r-1)*n1+1:(r-1)*n1+n1,(c-1)*n2+1:(c-1)*n2+n2) = zeros(n1,n2);
            end
            figscount = figscount+1;
        end
    end
