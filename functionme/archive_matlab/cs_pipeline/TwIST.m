function [x,x_debias,objective,times,debias_start,mses,max_svd] = ...
         TwIST(y,A,tau,varargin)
%
% Usage:
% [x,x_debias,objective,times,debias_start,mses] = TwIST(y,A,tau,varargin)
%
% This function solves the regularization problem 
%
%     arg min_x = 0.5*|| y - A x ||_2^2 + tau phi( x ), 
%
% where A is a generic matrix and phi(.) is a regularizarion 
% function  such that the solution of the denoising problem 
%
%     Psi_tau(y) = arg min_x = 0.5*|| y - x ||_2^2 + tau \phi( x ), 
%
% is known. 
% 
% For further details about the TwIST algorithm, see the paper:
%
% J. Bioucas-Dias and M. Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image 
% Restoration",  IEEE Transactions on Image processing, 2007.
% 
% and
% 
% J. Bioucas-Dias and M. Figueiredo, "A Monotonic Two-Step 
% Algorithm for Compressive Sensing and Other Ill-Posed 
% Inverse Problems", submitted, 2007.
%
% Authors: Jose Bioucas-Dias and Mario Figueiredo, October, 2007.
% 
% Please check for the latest version of the code and papers at
% www.lx.it.pt/~bioucas/TwIST
%
% -----------------------------------------------------------------------
% Copyright (2007): Jose Bioucas-Dias and Mario Figueiredo
% 
% TwIST is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% 
%  ===== Required inputs =============
%
%  y: 1D vector or 2D array (image) of observations
%     
%  A: if y and x are both 1D vectors, A can be a 
%     k*n (where k is the size of y and n the size of x)
%     matrix or a handle to a function that computes
%     products of the form A*v, for some vector v.
%     In any other case (if y and/or x are 2D arrays), 
%     A has to be passed as a handle to a function which computes 
%     products of the form A*x; another handle to a function 
%     AT which computes products of the form A'*x is also required 
%     in this case. The size of x is determined as the size
%     of the result of applying AT.
%
%  tau: regularization parameter, usually a non-negative real 
%       parameter of the objective  function (see above). 
%  
%
%  ===== Optional inputs =============
%  
%  'Psi' = denoising function handle; handle to denoising function
%          Default = soft threshold.
%
%  'Phi' = function handle to regularizer needed to compute the objective
%          function.
%          Default = ||x||_1
%
%  'lambda' = lam1 parameters of the  TwIST algorithm:
%             Optimal choice: lam1 = min eigenvalue of A'*A.
%             If min eigenvalue of A'*A == 0, or unknwon,  
%             set lam1 to a value much smaller than 1.
%
%             Rule of Thumb: 
%                 lam1=1e-4 for severyly ill-conditioned problems
%                 lam1=1e-2 for mildly  ill-conditioned problems
%                 lam1=1    for A unitary direct operators
%
%             Default: lam1 = 0.04.
%
%             Important Note: If (max eigenvalue of A'*A) > 1,
%             the algorithm may diverge. This is  be avoided 
%             by taking one of the follwoing  measures:
% 
%                1) Set 'Monontone' = 1 (default)
%                  
%                2) Solve the equivalenve minimization problem
%
%             min_x = 0.5*|| (y/c) - (A/c) x ||_2^2 + (tau/c^2) \phi( x ), 
%
%             where c > 0 ensures that  max eigenvalue of (A'A/c^2) <= 1.
%
%   'alpha' = parameter alpha of TwIST (see ex. (22) of the paper)         
%             Default alpha = alpha(lamN=1, lam1)
%   
%   'beta'  =  parameter beta of twist (see ex. (23) of the paper)
%              Default beta = beta(lamN=1, lam1)            
% 
%  'AT'    = function handle for the function that implements
%            the multiplication by the conjugate of A, when A
%            is a function handle. 
%            If A is an array, AT is ignored.
%
%  'StopCriterion' = type of stopping criterion to use
%                    0 = algorithm stops when the relative 
%                        change in the number of non-zero 
%                        components of the estimate falls 
%                        below 'ToleranceA'
%                    1 = stop when the relative 
%                        change in the objective function 
%                        falls below 'ToleranceA'
%                    2 = stop when the relative norm of the difference between 
%                        two consecutive estimates falls below toleranceA
%                    3 = stop when the objective function 
%                        becomes equal or less than toleranceA.
%                    Default = 1.
%
%  'ToleranceA' = stopping threshold; Default = 0.01
% 
%  'Debias'     = debiasing option: 1 = yes, 0 = no.
%                 Default = 0.
%                 
%                 Note: Debiasing is an operation aimed at the 
%                 computing the solution of the LS problem 
%
%                         arg min_x = 0.5*|| y - A' x' ||_2^2 
%
%                 where A' is the  submatrix of A obatained by
%                 deleting the columns of A corresponding of components
%                 of x set to zero by the TwIST algorithm
%                 
%
%  'ToleranceD' = stopping threshold for the debiasing phase:
%                 Default = 0.0001.
%                 If no debiasing takes place, this parameter,
%                 if present, is ignored.
%
%  'MaxiterA' = maximum number of iterations allowed in the
%               main phase of the algorithm.
%               Default = 1000
%
%  'MiniterA' = minimum number of iterations performed in the
%               main phase of the algorithm.
%               Default = 5
%
%  'MaxiterD' = maximum number of iterations allowed in the
%               debising phase of the algorithm.
%               Default = 200
%
%  'MiniterD' = minimum number of iterations to perform in the
%               debiasing phase of the algorithm.
%               Default = 5
%
%  'Initialization' must be one of {0,1,2,array}
%               0 -> Initialization at zero. 
%               1 -> Random initialization.
%               2 -> initialization with A'*y.
%               array -> initialization provided by the user.
%               Default = 0;
%
%  'Monotone' = enforce monotonic decrease in f. 
%               any nonzero -> enforce monotonicity
%               0 -> don't enforce monotonicity.
%               Default = 1;
%
%  'Sparse'   = {0,1} accelarates the convergence rate when the regularizer 
%               Phi(x) is sparse inducing, such as ||x||_1.
%               Default = 1
%               
%             
%  'True_x' = if the true underlying x is passed in 
%                this argument, MSE evolution is computed
%
%
%  'Verbose'  = work silently (0) or verbosely (1)
%
% ===================================================  
% ============ Outputs ==============================
%   x = solution of the main algorithm
%
%   x_debias = solution after the debiasing phase;
%                  if no debiasing phase took place, this
%                  variable is empty, x_debias = [].
%
%   objective = sequence of values of the objective function
%
%   times = CPU time after each iteration
%
%   debias_start = iteration number at which the debiasing 
%                  phase started. If no debiasing took place,
%                  this variable is returned as zero.
%
%   mses = sequence of MSE values, with respect to True_x,
%          if it was given; if it was not given, mses is empty,
%          mses = [].
%
%   max_svd = inverse of the scaling factor, determined by TwIST,
%             applied to the direct operator (A/max_svd) such that
%             every IST step is increasing.
% ========================================================

% 【在CS链路中的位置】
%   TwIST是压缩感知重建的核心优化求解器，被AS_reconstructCS调用。
%   求解问题：min 0.5||y - A(x)||² + tau·Phi(x)
%   其中 y = g_vec（观测全息图实向量），A = MyForwardOperatorPropagation
%   Phi = MyTVphi（TV正则化），Psi = MyTVpsi（TV去噪近端算子）
%   调用链：AS_reconstructCS → TwIST(g_vec, A, tau, 'AT', AT, 'Psi', Psi, 'Phi', Phi, ...)
%
% 【算法原理】
%   TwIST (Two-step Iterative Shrinkage/Thresholding) 是一种两步迭代算法，
%   相比单步IST具有更快的收敛速度。核心迭代：
%     1. 计算梯度：grad = AT(y - A(x))
%     2. IST步：x_ist = Psi(x + grad/max_svd, tau/max_svd)
%     3. 两步外推：x_new = (alpha-beta)*x_prev + (1-alpha)*x_prev2 + beta*x_ist
%     4. 若单调性被破坏则回退到IST步
%
% 【输入】
%   y     - 观测数据实向量 [2·Nx·Ny × 1]，全息图的[real;imag]编码
%   A     - 正向算子函数句柄 @(f) → MyForwardOperatorPropagation
%   tau   - 正则化参数（标量）
%   可选参数通过varargin传入（键值对形式）
%
% 【输出】
%   x           - 解向量 [2·Nx·Ny·Nz × 1]，重建散射势的[real;imag]编码
%   x_debias    - 去偏后的解（若Debias=1）
%   objective   - 目标函数值序列
%   times       - 每次迭代的CPU时间
%   debias_start - 去偏阶段开始的迭代号
%   mses        - MSE序列（若提供True_x）
%   max_svd     - A'A最大特征值的估计（用于步长缩放）

%--------------------------------------------------------------
% 检查必需参数数量（y, A, tau 共3个）
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 3
     error('Wrong number of required parameters');
end
%--------------------------------------------------------------
% 设置可选参数的默认值
%--------------------------------------------------------------
stopCriterion = 1;      % 停止准则：1=目标函数相对变化
tolA = 0.01;            % 主循环收敛容差
debias = 0;             % 是否执行去偏阶段
maxiter = 1000;         % 主循环最大迭代次数
maxiter_debias = 200;   % 去偏阶段最大迭代次数
miniter = 5;            % 主循环最小迭代次数
miniter_debias = 5;     % 去偏阶段最小迭代次数
init = 0;               % 初始化方式：0=零，1=随机，2=AT(y)
enforceMonotone = 1;    % 是否强制目标函数单调下降
compute_mse = 0;        % 是否计算MSE（需要True_x）
plot_ISNR = 0;          % 是否绘制ISNR
AT = 0;                 % 伴随算子句柄（默认0，后续检查）
verbose = 1;            % 是否打印迭代信息
alpha = 0;              % TwIST参数alpha（0表示自动计算）
beta  = 0;              % TwIST参数beta（0表示自动计算）
sparse = 1;             % 是否启用稀疏加速
tolD = 0.001;           % 去偏阶段收敛容差
phi_l1 = 0;             % 标记是否使用L1正则（默认）
psi_ok = 0;             % 标记Psi函数句柄是否验证通过
% 默认特征值范围：lam1=A'A最小特征值估计，lamN=A'A最大特征值估计（归一化为1）
lam1=1e-4;   lamN=1;
% 

% 常量与内部变量
for_ever = 1;           % 用于内部while循环的无限循环标志
% max_svd: 算子A的最大奇异值的上界估计，用于缩放步长
max_svd = 1;

% 初始化可能不被计算的输出变量
debias_start = 0;
x_debias = [];
mses = [];

%--------------------------------------------------------------
% 读取可选参数（键值对形式）
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch upper(varargin{i})
     case 'LAMBDA'
       lam1 = varargin{i+1};
     case 'ALPHA'
       alpha = varargin{i+1};
     case 'BETA'
       beta = varargin{i+1};
     case 'PSI'
       psi_function = varargin{i+1};   % TV去噪近端算子句柄
     case 'PHI'
       phi_function = varargin{i+1};   % TV正则化函数值句柄
     case 'STOPCRITERION'
       stopCriterion = varargin{i+1};
     case 'TOLERANCEA'       
       tolA = varargin{i+1};
     case 'TOLERANCED'
       tolD = varargin{i+1};
     case 'DEBIAS'
       debias = varargin{i+1};
     case 'MAXITERA'
       maxiter = varargin{i+1};
     case 'MAXIRERD'
       maxiter_debias = varargin{i+1};
     case 'MINITERA'
       miniter = varargin{i+1};
     case 'MINITERD'
       miniter_debias = varargin{i+1};
     case 'INITIALIZATION'
       if prod(size(varargin{i+1})) > 1   % 用户提供了初始x
	 init = 33333;    % 特殊标志，表示用户提供了初始值
	 x = varargin{i+1};
       else 
	 init = varargin{i+1};
       end
     case 'MONOTONE'
       enforceMonotone = varargin{i+1};
     case 'SPARSE'
       sparse = varargin{i+1};
     case 'TRUE_X'
       compute_mse = 1;
       true = varargin{i+1};
        size(true)
        size(y)
       if prod(double((size(true) == size(y))))
           plot_ISNR = 1;
       end
     case 'AT'
       AT = varargin{i+1};   % 伴随算子函数句柄
     case 'VERBOSE'
       verbose = varargin{i+1};
     otherwise
      % 未识别的参数名
      error(['Unrecognized option: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%


% ====== 计算TwIST算法参数alpha和beta ======
% rho0 与 A'A 的条件数相关，用于确定最优的两步外推参数
rho0 = (1-lam1/lamN)/(1+lam1/lamN);
if alpha == 0 
    % alpha的最优值（论文公式(22)），当lamN=1时简化
    alpha = 2/(1+sqrt(1-rho0^2));
end
if  beta == 0 
    % beta的最优值（论文公式(23)）
    beta  = alpha*2/(lam1+lamN);
end


% 验证停止准则合法性
if (sum(stopCriterion == [0 1 2 3])==0)
   error(['Unknwon stopping criterion']);
end

% 如果A是函数句柄，必须同时提供AT（伴随算子）
if isa(A, 'function_handle') & ~isa(AT,'function_handle')
   error(['The function handle for transpose of A is missing']);
end 


% 如果A是矩阵（非函数句柄），则将其转换为函数句柄形式
% 以便后续代码统一使用函数句柄接口
if ~isa(A, 'function_handle')
   AT = @(x) reshape(A'*x(:),[64 64]);
   A = @(x) reshape(A*x(:),[size(y,1) size(y,2)]);
end
% 从此处开始，A和AT始终为函数句柄

% 预计算 AT(y)，后续多处使用（伴随反投影）
Aty = AT(y);
% psi_function(Aty,tau)

% ====== 验证Psi去噪函数句柄 ======
% 如果提供了Psi，验证其是否为合法函数句柄且能正确执行
if exist('psi_function','var')
   if isa(psi_function,'function_handle')
       try  % 用Aty测试Psi是否能正常工作
            dummy = psi_function(Aty,tau); 
            psi_ok = 1;
      catch ME
         error('TwIST:InvalidPsi', ...
               'Something is wrong with function handle for psi: %s', ...
               ME.message)
      end
   else
      error(['Psi does not seem to be a valid function handle']);
   end
else %if nothing was given, use soft thresholding
   % 未提供Psi时，默认使用软阈值算子
   psi_function = @(x,tau) soft(x,tau);
end

% ====== 验证Phi正则化函数句柄 ======
% 如果Psi存在，Phi也必须存在（用于计算目标函数值）
if (psi_ok == 1)
   if exist('phi_function','var')
      if isa(phi_function,'function_handle')
         try  % 测试Phi是否能正常工作
              dummy = phi_function(Aty); 
         catch ME
           error('TwIST:InvalidPhi', ...
                 'Something is wrong with function handle for phi: %s', ...
                 ME.message)
         end
      else
        error(['Phi does not seem to be a valid function handle']);
      end
   else
      error(['If you give Psi you must also give Phi']); 
   end
else  % 若未提供Psi和Phi，默认使用L1范数作为正则化
   phi_function = @(x) sum(abs(x(:))); 
   phi_l1 = 1;  % 标记正在使用L1正则
end
    

%--------------------------------------------------------------
% 初始化解向量x
%--------------------------------------------------------------
switch init
    case 0   % 初始化为零向量，通过AT确定x的尺寸
       x = AT(zeros(size(y), 'like', y));
    case 1   % 随机初始化
       x = randn(size(AT(zeros(size(y), 'like', y))), 'like', y);
    case 2   % 用AT(y)初始化（伴随反投影），通常是最实用的初始化
       x = Aty; 
    case 33333
       % 用户提供了初始x，检查其尺寸与A是否兼容
       if size(A(x)) ~= size(y)
          error(['Size of initial x is not compatible with A']); 
       end
    otherwise
       error(['Unknown ''Initialization'' option']);
end

% 检查tau的尺寸：标量或与x同尺寸
if prod(size(tau)) > 1
   try,
      dummy = x.*tau;
   catch,
      error(['Parameter tau has wrong dimensions; it should be scalar or size(x)']),
   end
end
      
% 如果提供了真实解true_x，检查其尺寸
if compute_mse & (size(true) ~= size(x))  
   error(['Initial x has incompatible size']); 
end


% ====== L1正则化下的零解检测 ======
% 如果tau足够大（>= max|AT(y)|），软阈值会直接将所有分量收缩为零
% 此时最优解就是零向量，无需迭代
if phi_l1
   max_tau = toCpuScalar(max(abs(Aty(:))));
   if (tau >= max_tau) && (psi_ok==0)
      x = zeros(size(Aty), 'like', Aty);
      objective(1) = toCpuScalar(0.5*(y(:)'*y(:)));
      times(1) = 0;
      if compute_mse
        mses(1) = toCpuScalar(sum(true(:).^2));
      end
      return
   end
end


% 统计初始解中非零元素的数量
nz_x = (x ~= 0.0);
num_nz_x = toCpuScalar(sum(nz_x(:)));

% 计算并存储初始目标函数值：0.5||y-Ax||² + tau·Phi(x)
resid =  y-A(x);
prev_f = toCpuScalar(0.5*(resid(:)'*resid(:)) + tau*phi_function(x));


% 启动计时
t0 = cputime;

times(1) = cputime - t0;
objective(1) = prev_f;

if compute_mse
   mses(1) = toCpuScalar(sum(sum((x-true).^2)));
end

cont_outer = 1;   % 外层循环控制标志
iter = 1;

if verbose
    fprintf(1,'\nInitial objective = %10.6e,  nonzeros=%7d\n',...
        prev_f,num_nz_x);
end

% 控制IST和TwIST步骤计数的变量
IST_iters = 0;    % 连续IST步计数
TwIST_iters = 0;  % TwIST步计数

% 初始化历史解（两步外推需要前两步的解）
xm2=x;  % x_{k-2}
xm1=x;  % x_{k-1}

%--------------------------------------------------------------
% TwIST主迭代循环
%--------------------------------------------------------------
while cont_outer
    % 计算梯度：grad = AT(resid) = AT(y - A(x))
    grad = AT(resid);
    while for_ever
        % ====== IST步：带近端算子的梯度下降 ======
        % x_ist = Psi(xm1 + grad/max_svd, tau/max_svd)
        % max_svd用于缩放梯度步长，确保收敛
        x = psi_function(xm1 + grad/max_svd,tau/max_svd);
        if (IST_iters >= 2) | ( TwIST_iters ~= 0)
            % ====== 两步外推（TwIST核心步骤）======
            % 当已有足够IST步数或之前已进入TwIST模式时，执行两步外推
            if sparse
                % 稀疏加速：将当前为零的位置在过去解中也置零
                % 防止已被收缩为零的分量通过外推被重新激活
                mask = (x ~= 0);
                xm1 = xm1.* mask;
                xm2 = xm2.* mask;
            end
            % 两步外推公式：x_new = (alpha-beta)*x_{k-1} + (1-alpha)*x_{k-2} + beta*x_ist
            xm2 = (alpha-beta)*xm1 + (1-alpha)*xm2 + beta*x;
            % 计算外推后的残差和目标函数
            resid = y-A(xm2);
            f = toCpuScalar(0.5*(resid(:)'*resid(:)) + tau*phi_function(xm2));
            if (f > prev_f) & (enforceMonotone)
                % 单调性被破坏：回退到IST步
                TwIST_iters = 0;  % do a IST iteration if monotonocity fails
            else
                % 单调性满足：接受TwIST步
                TwIST_iters = TwIST_iters+1; % TwIST iterations
                IST_iters = 0;
                x = xm2;
                % 定期衰减max_svd以加速收敛（减少保守的步长缩放）
                if mod(TwIST_iters,10000) == 0
                    max_svd = 0.9*max_svd;
                end
                break;  % 退出内层while循环
            end
        else
            % ====== 纯IST步（前两次迭代或TwIST回退后）======
            resid = y-A(x);
            f = toCpuScalar(0.5*(resid(:)'*resid(:)) + tau*phi_function(x));
            if f > prev_f
                % IST步单调性也失败，说明max_svd估计不足
                % 增大max_svd（减小步长），防止A'A的最大特征值>1导致发散
                max_svd = 2*max_svd;
                if verbose
                    fprintf('Incrementing S=%2.2e\n',max_svd)
                end
                IST_iters = 0;
                TwIST_iters = 0;
            else
                % IST步成功
                TwIST_iters = TwIST_iters + 1;
                break;  % 退出内层while循环
            end
        end
    end

    % 更新历史解
    xm2 = xm1;  % x_{k-2} ← 旧的x_{k-1}
    xm1 = x;     % x_{k-1} ← 当前的x


    % ====== 更新非零分量统计 ======
    nz_x_prev = nz_x;
    nz_x = (x~=0.0);
    num_nz_x = toCpuScalar(sum(nz_x(:)));
    num_changes_active = toCpuScalar(sum(nz_x(:)~=nz_x_prev(:)));

    % ====== 计算停止准则 ======
    % 至少执行miniter次，至多执行maxiter次
    switch stopCriterion
        case 0,
            % 基于非零分量数变化
            criterion =  num_changes_active;
        case 1,
            % 基于目标函数相对变化
            criterion = abs(f-prev_f)/prev_f;
        case 2,
            % 基于解的相对范数变化
            criterion = (norm(x(:)-xm1(:))/norm(x(:)));
        case 3,
            % 基于目标函数绝对值
            criterion = f;
        otherwise,
            error(['Unknwon stopping criterion']);
    end
    criterion = toCpuScalar(criterion);
    cont_outer = ((iter <= maxiter) && (criterion > tolA));  % 未达最大迭代且未收敛
    if iter <= miniter
        cont_outer = 1;  % 强制至少执行miniter次
    end



    iter = iter + 1;
    prev_f = f;
    objective(iter) = f;          % 记录目标函数值
    times(iter) = cputime-t0;     % 记录累计CPU时间

    if compute_mse
        err = true - x;
        mses(iter) = toCpuScalar(err(:)'*err(:));
    end

    % 打印迭代信息
    if verbose
        if plot_ISNR
            fprintf(1,'Iteration=%4d, ISNR=%4.5e  objective=%9.5e, nz=%7d, criterion=%7.3e\n',...
                iter, 10*log10(sum((y(:)-true(:)).^2)/sum((x(:)-true(:)).^2) ), ...
                f, num_nz_x, criterion/tolA);
        else
            fprintf(1,'Iteration=%4d, objective=%9.5e, nz=%7d,  criterion=%7.3e\n',...
                iter, f, num_nz_x, criterion/tolA);
        end
    end

%     figure(999);imagesc(plotdatacube(x));colormap gray;axis image;colorbar;drawnow;
end
%--------------------------------------------------------------
% 主循环结束
%--------------------------------------------------------------

% 打印最终结果摘要
if verbose
    fprintf(1,'\nFinished the main algorithm!\nResults:\n')
    fprintf(1,'||A x - y ||_2 = %10.3e\n',toCpuScalar(resid(:)'*resid(:)))
    fprintf(1,'||x||_1 = %10.3e\n',toCpuScalar(sum(abs(x(:)))))
    fprintf(1,'Objective function = %10.3e\n',f);
    fprintf(1,'Number of non-zero components = %d\n',num_nz_x);
    fprintf(1,'CPU time so far = %10.3e\n', times(iter));
    fprintf(1,'\n');
end


%--------------------------------------------------------------
% 去偏阶段（Debiasing）
% 当Debias=1时，在TwIST主循环找到稀疏支撑后，
% 用共轭梯度法(CG)在非零支撑上求解无正则化的最小二乘问题，
% 以消除L1/TV正则化引入的幅值偏差。
%--------------------------------------------------------------
if debias
    if verbose
        fprintf(1,'\n')
        fprintf(1,'Starting the debiasing phase...\n\n')
    end

    x_debias = x;
    zeroind = (x_debias~=0);  % 非零元素的掩模（支撑集）
    cont_debias_cg = 1;
    debias_start = iter;

    % 计算初始残差
    resid = A(x_debias);
    resid = resid-y;
    resid_prev = eps*ones(size(resid));

    rvec = AT(resid);

    % 用支撑集掩模屏蔽零位置（固定零分量为零）
    rvec = rvec .* zeroind;
    rTr_cg = rvec(:)'*rvec(:);

    % 设置去偏阶段的残差收敛阈值
    tol_debias = tolD * (rvec(:)'*rvec(:));

    % 初始化CG搜索方向
    pvec = -rvec;

    % CG主循环
    while cont_debias_cg

        % 计算 A*p = AT(A(pvec))，只在支撑集上
        RWpvec = A(pvec);
        Apvec = AT(RWpvec);

        % 掩模屏蔽
        Apvec = Apvec .* zeroind;

        % CG步长alpha
        alpha_cg = rTr_cg / (pvec(:)'* Apvec(:));

        % 更新解、残差和梯度
        x_debias = x_debias + alpha_cg * pvec;
        resid = resid + alpha_cg * RWpvec;
        rvec  = rvec  + alpha_cg * Apvec;

        % 计算新的rTr，更新CG搜索方向
        rTr_cg_plus = rvec(:)'*rvec(:);
        beta_cg = rTr_cg_plus / rTr_cg;
        pvec = -rvec + beta_cg * pvec;

        rTr_cg = rTr_cg_plus;

        iter = iter+1;

        objective(iter) = 0.5*(resid(:)'*resid(:)) + ...
            tau*phi_function(x_debias(:));
        times(iter) = cputime - t0;

        if compute_mse
            err = true - x_debias;
            mses(iter) = (err(:)'*err(:));
        end

        % CG收敛判断：基于残差准则
        if verbose
            fprintf(1,' Iter = %5d, debias resid = %13.8e, convergence = %8.3e\n', ...
                iter, resid(:)'*resid(:), rTr_cg / tol_debias);
        end
        cont_debias_cg = ...
            (iter-debias_start <= miniter_debias )| ...
            ((rTr_cg > tol_debias) & ...
            (iter-debias_start <= maxiter_debias));

    end
    if verbose
        fprintf(1,'\nFinished the debiasing phase!\nResults:\n')
        fprintf(1,'||A x - y ||_2 = %10.3e\n',toCpuScalar(resid(:)'*resid(:)))
        fprintf(1,'||x||_1 = %10.3e\n',toCpuScalar(sum(abs(x(:)))))
        fprintf(1,'Objective function = %10.3e\n',f);
        nz = (x_debias~=0.0);
        fprintf(1,'Number of non-zero components = %d\n',toCpuScalar(sum(nz(:))));
        fprintf(1,'CPU time so far = %10.3e\n', times(iter));
        fprintf(1,'\n');
    end
end

if compute_mse
   mses = mses/length(true(:));
end


%--------------------------------------------------------------
% soft  软阈值算子，适用于实数和复数
% 当未提供自定义Psi时，作为默认的近端算子
%--------------------------------------------------------------
function y = soft(x,T)
%y = sign(x).*max(abs(x)-tau,0);
y = max(abs(x) - T, 0);      % 幅值收缩：(|x|-T)+ 
y = y./(y+T) .* x;           % 保持相位，幅值按 (|x|-T)/|x| 缩放

function s = toCpuScalar(v)
% toCpuScalar  将GPU标量回收到CPU并转为double
%   确保控制流判断中使用的是CPU上的实数标量，避免GPU数据类型问题
if isa(v, 'gpuArray')
    v = gather(v);    % GPU → CPU
end
s = double(v);
if ~isreal(s)
    s = real(s);      % 取实部，避免复数残留
end
if numel(s) ~= 1
    error('TwIST:ExpectedScalar', 'Expected scalar value for control flow.');
end
