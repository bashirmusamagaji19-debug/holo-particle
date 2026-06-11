%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
数值范围映射 (Value Range Mapping to [min_val, max_val])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
【用途】
  将图像像素值线性映射到指定范围 [min_val, max_val]。
  支持多通道（如RGB）图像，各通道使用相同的全局最大最小值进行映射，
  保持通道间的相对比例关系。

【数学公式】
  令 maxI = max(oriImage(:))，minI = min(oriImage(:))
  线性映射系数：
    coeA = (max_val - min_val) / (maxI - minI)
    coeB = min_val - minI · coeA
  输出：outImage = oriImage · coeA + coeB

  特殊情况：当 maxI == minI 时（如全黑图像），直接缩放：
    outImage = oriImage · (max_val / minI)

【输入】
  oriImage - 输入图像 (Ny × Nx × Chan)，可为灰度或彩色
  min_val  - 映射范围下界
  max_val  - 映射范围上界

【输出】
  outImage - 映射后图像 (Ny × Nx × Chan)，double 类型

【作者】Ni Chen (ni_chen@163.com), May 2012
【版权】GNU General Public License v2+
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
function outImage = Toxy(oriImage, min_val, max_val)

    % 获取图像尺寸：行数、列数、通道数
    [Ny, Nx, Chan] = size(oriImage);

    % 预分配输出矩阵
    outImage = zeros(Ny, Nx, Chan);

    % 计算输入图像的全局最大值和最小值（跨所有通道）
    maxI = max(max(double(oriImage(:))));
    minI = min(min(double(oriImage(:))));
    
    % 逐通道进行映射
    for c = 1:Chan
        tempImg = oriImage(:,:,c);

        % 特殊情况：所有像素值相同
        if maxI == minI
            % 若最小值非零，直接按比例缩放到 max_val
            if minI~=0
                outImage(:,:,c) = double(tempImg)*max_val/minI;
            end
            % 若最小值为零（全黑图像），输出保持全零
        else
            % 一般情况：线性映射 y = coeA·x + coeB
            coeA = (max_val-min_val)/double(maxI-minI);
            coeB = min_val-minI*coeA;
            outImage(:,:,c) = double(tempImg).*coeA+coeB;
        end
    end 
end
