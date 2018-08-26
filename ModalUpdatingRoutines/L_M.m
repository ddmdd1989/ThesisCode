function [x, RESNORM, RESIDUAL, OUTPUT] = L_M(rfun, x0, options)
% function [x, RESNORM, RESIDUAL, OUTPUT] = L_M(rFun, x0, options)
%
%	(c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function uses Levenberg-Marquardt algorithm to solve a box-constrained 
% optimization problem. The objective is to minimize the 2-norm of an m x 1 
% residual vector for an n x 1 optimization variable. 
%
% Inputs:
%   [y,Jac] = rFun(x) 
%       y: m x 1 objective residual vector r. The vector contains the m number
%           of residuals whose sum of squares are to be minimized.
%       Jac: The m x n Jacobian matrix of the residuals with respect to
%           the optimzation variables (dr_i / dx_j)
%   x0: n x 1 initial starting point of the optimization variables
%	options: user defined options for Gauss-Newton optimization process
%		maxNumItr: Maximum number of iterations for the Levenberg-Marquardt
%		optimization; (default: 1e4)
%   	tolGrad: Termination tolerance on the gradient at current values
%				 of optimization variables; (default: 1e-3)
%   	tolX: Termination tolerance on the value change of the optimization
%      		  variables between two iterations; (default: 1e-5)
%   	tolFun: Termination tolerance on the value change of the objective
%      			function between two iterations;(default: 1e-5)
% Outputs:
% 	x: optimal values of optimization variables
% 	RESNORM: 2-norm of the optimal residual vector
% 	RESIDUAL: optimal residual vector
%   OUTPUT: additional information of the optimization solution
%       grad: gradient of the sum of squares of r_i over x
%       numIter: Number of iterations of Gauss-Newton optmization process.
%       exitflag: exit conditions with following four values
%    		1: TolGrad termination criterion
%       	2: TolX termination criterion
%           3: TolFun termination criterion
%    		4: Exceed the specified maximum number of iterations (maxNumItr)

%% Options default values
tolGrad = 1e-3;
tolX = 1e-5;
tolFun = 1e-5;
maxNumItr = 1e4;
lb = -Inf * ones(length(x0),1);
ub = Inf*ones(length(x0),1);

%% Use user-defined options if exists
if(nargin == 3)
    if(isfield(options,'tolGrad'))
        tolGrad = options.tolGrad;
    end
    if(isfield(options,'tolX'))
        tolX = options.tolX;
    end
    
    if(isfield(options,'tolFun'))
        tolFun = options.tolFun;
    end
    if(isfield(options,'maxNumItr'))
        maxNumItr = options.maxNumItr;
    end
    if(isfield(options,'lb'))
        lb = options.lb;
    end
    if(isfield(options,'ub'))
        ub = options.ub;
    end
   
end



header = sprintf(['\n                                  Norm of      First-order \n',...
    ' Iteration          f(x)          step          optimality          Lambda\n']);

fprintf(header);

formatstr = ' %5.0f     %13.6g  %13.6g   %12.3g      %12.3g\n';

alpha = 0.01;
lambda = 0.01;
x = x0;
bk_tr = 0;
for k = 1 : maxNumItr
    initLambda = lambda;
    [obj, A] = feval(rfun, x);
    if (norm(2 * A' * obj) < tolGrad)
        newobj = obj;
        EXITFLAG = 1;
        break;
    end
    
    AtA = A' * A;
    Atr = A' * obj;
    
    reg_term = eye(size(AtA,1));
    
    v = -(AtA + lambda * reg_term) \ Atr;
   
    temp1 =(x + v < lb);
    temp2 =(x + v > ub);
    % Ensure the updated value is within the bound
    while (~isempty(find(temp1,1)) || ~isempty(find(temp2,1)))
        lambda = 10 * lambda;
        v = -(AtA + lambda * reg_term) \ Atr;
        temp1 =(x + v < lb);
        temp2 =(x + v > ub);
        bk_tr = bk_tr + 1;
    end
    if(~isempty(find(v == Inf, 1)) || ~isempty(find(isnan(v), 1)))
        EXITFLAG = 1;
        break;
    end
        
    newx = x + v;
    
    newobj = feval(rfun,newx);
    
    
    while (norm(newobj)^2 > norm(obj)^2 + alpha * (2*obj'*A*v)  * initLambda /lambda)
        lambda = 10 * lambda;
        v = -(AtA + lambda * reg_term) \ Atr;
        newx = x + v;
        newobj = feval(rfun,newx);
        bk_tr = bk_tr + 1;
    end
        
    if(bk_tr == 0) 
        lambda = 0.1 * lambda;
    end
    bk_tr = 0;
    

    x = newx;
    
    if(norm(v) < tolX)
        EXITFLAG = 2;
        break;
    end
    if(abs(norm(newobj)^2 - norm(obj)^2) < tolFun)
        EXITFLAG = 3;
        break;
    end
    fprintf(formatstr,k,norm(newobj)^2,norm(v),norm(2*A'*obj),lambda)
end

if(k == maxNumItr)
    EXITFLAG = 4;
end

RESNORM = norm(newobj)^2;
RESIDUAL = newobj;
OUTPUT.gradient = 2 * A' * newobj;
OUTPUT.iterations = k;
OUTPUT.exitflag = EXITFLAG;

end

