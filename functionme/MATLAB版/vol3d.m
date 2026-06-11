function [model] = vol3d(varargin)
%VOL3D 底层三维体绘制函数（基于正交切片纹理映射）
%
%   【用途】
%   使用正交平面2D纹理映射技术，在OpenGL中实现3D数据的体绘制渲染。
%   通过在x/y/z三个方向堆叠半透明纹理切片来模拟体积渲染效果。
%   需要较快的OpenGL硬件支持以获得流畅体验。
%
%   【在工程中的位置】
%   这是3D渲染的底层函数。show3d.m 和 GUI 中的3D显示均依赖此函数。
%   属于基础设施层，被上层显示函数调用。
%
%   【调用方式】
%   H = vol3d('CData',data)              从3D数据创建体渲染对象，返回结构体句柄
%   vol3d(...,'Parent',axH)              指定父坐标轴
%   vol3d(...,'texture','2D')            默认模式：仅渲染最接近视角的正交切片，需手动刷新
%   vol3d(...,'texture','3D')            同时渲染x/y/z三个方向切片，无需刷新但需更强GPU
%   vol3d(H)                             刷新视图（'2D'模式下旋转后需调用）
%
%   【注意事项】
%   - 可配合 vol3dtool 编辑色图和透明度图
%   - 可用 interp3 对输入数据进行插值以调整渲染分辨率
%   - 本项目统一使用 '3D' 纹理模式

% Copyright Joe Conti, 2004

%% ---- 解析输入参数 ----
if isstruct(varargin{1})
    % 传入已有model结构体时（用于刷新视图）
    model = varargin{1};
    if length(varargin) > 1
       varargin = {varargin{2:end}};
    end
else
    % 首次调用，使用默认模型
    model = localGetDefaultModel;
end


if length(varargin)>1
  % 解析名值对参数（CData / Parent / texture）
  for n = 1:2:length(varargin)
    switch(lower(varargin{n}))
        case 'cdata'
            model.cdata = varargin{n+1};     % 体数据
        case 'parent'
            model.parent = varargin{n+1};    % 父坐标轴句柄
        case 'texture'
            model.texture = varargin{n+1};   % 纹理模式：'2D' 或 '3D'
    end
    
  end
end

if isempty(model.parent)
    model.parent = gca;  % 默认使用当前坐标轴
end

% 设置坐标轴为3D视角模式
ax = model.parent;
axis(ax,'vis3d');
axis(ax,'tight');

% 执行绘制
[model] = local_draw(model);


%------------------------------------------%
% 获取默认模型参数
%------------------------------------------%
function [model] = localGetDefaultModel

model.cdata = [];      % 体数据（待绘制）
model.xdata = [];      % x轴数据范围
model.ydata = [];      % y轴数据范围
model.zdata = [];      % z轴数据范围
model.parent = [];     % 父坐标轴句柄
model.handles = [];    % 图形对象句柄数组
model.texture = '2D';  % 默认2D纹理模式

%------------------------------------------%
% 核心绘制函数：根据相机方向和纹理模式创建正交切片
%------------------------------------------%
function [model,ax] = local_draw(model)

cdata = model.cdata; 
siz = size(cdata);  % siz = [行数, 列数, 层数] = [ny, nx, nz]

% 若未指定坐标范围，则使用数据索引作为默认范围
if isempty(model.xdata)
    model.xdata = [0 siz(2)];  % x方向：0 到 nx
end
if isempty(model.ydata)
    model.ydata = [0 siz(1)];  % y方向：0 到 ny
end
if isempty(model.zdata)
    model.zdata = [0 siz(3)];  % z方向：0 到 nz
end

% 清除之前的图形对象
try,
   delete(model.handles);
end

ax = model.parent;
% 获取当前相机方向，确定最接近的正交平面
cam_dir = camtarget(ax) - campos(ax);  % 相机方向向量
[m,ind] = max(abs(cam_dir));           % ind=1/2/3 分别对应最接近x/y/z方向

% 清除坐标轴中已有的vol3d对象
h = findobj(ax,'type','surface','tag','vol3d');
for n = 1:length(h)
  try,
     delete(h(n));
  end
end

% 判断是否使用3D纹理模式
is3DTexture = strcmpi(model.texture,'3D');
handle_ind = 1;
% 控制drawnow频率：减少drawnow调用次数以提高渲染吞吐量，同时保持UI响应
drawStride = max(1, floor(max(siz) / 12));

%% ---- 创建Z方向切片（沿z轴堆叠xy平面） ----
% 当相机最接近z方向(ind==3)或使用3D纹理模式时绘制
if(ind==3 || is3DTexture )    
  % 初始切片的四角坐标
  x = [model.xdata(1), model.xdata(2); model.xdata(1), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(1); model.zdata(1), model.zdata(1)];
  diff = model.zdata(2)-model.zdata(1);  % z方向总跨度
  delta = diff/size(cdata,3);            % 每层z切片间距
  for n = 1:size(cdata,3)

   slice = squeeze(cdata(:,:,n));       % 取第n层xy平面数据
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   % 设置纹理映射：颜色和透明度均由数据值驱动
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',slice,'facealpha','texturemap','tag','vol3d');
   z = z + delta;  % z坐标递增一层
   handle_ind = handle_ind + 1;
   % 每隔drawStride层刷新一次显示，平衡性能与响应性
   if mod(n, drawStride) == 0
      drawnow;
   end
  end

end

%% ---- 创建X方向切片（沿x轴堆叠yz平面） ----
% 当相机最接近x方向(ind==1)或使用3D纹理模式时绘制
if (ind==1 || is3DTexture ) 
  x = [model.xdata(1), model.xdata(1); model.xdata(1), model.xdata(1)];
  y = [model.ydata(1), model.ydata(1); model.ydata(2), model.ydata(2)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.xdata(2)-model.xdata(1);
  delta = diff/size(cdata,2);
  for n = 1:size(cdata,2)

   slice = squeeze(cdata(:,n,:));       % 取第n个yz平面数据
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',slice,'facealpha','texturemap','tag','vol3d');
   x = x + delta;  % x坐标递增一层
   handle_ind = handle_ind + 1;
   if mod(n, drawStride) == 0
      drawnow;
   end
  end
end

  
%% ---- 创建Y方向切片（沿y轴堆叠xz平面） ----
% 当相机最接近y方向(ind==2)或使用3D纹理模式时绘制
if (ind==2 || is3DTexture)
  x = [model.xdata(1), model.xdata(1); model.xdata(2), model.xdata(2)];
  y = [model.ydata(1), model.ydata(1); model.ydata(1), model.ydata(1)];
  z = [model.zdata(1), model.zdata(2); model.zdata(1), model.zdata(2)];
  diff = model.ydata(2)-model.ydata(1);
  delta = diff/size(cdata,1);
  for n = 1:size(cdata,1)

   slice = squeeze(cdata(n,:,:));       % 取第n个xz平面数据
   h(handle_ind) = surface(x,y,z,'Parent',ax);
   set(h(handle_ind),'cdatamapping','scaled','facecolor','texture','cdata',slice,...
	 'edgealpha',0,'alphadata',slice,'facealpha','texturemap','tag','vol3d');
   y = y + delta;  % y坐标递增一层
   handle_ind = handle_ind + 1;
   if mod(n, drawStride) == 0
      drawnow;
   end
  end
end

model.handles = h;  % 保存所有图形句柄，便于后续删除/刷新
