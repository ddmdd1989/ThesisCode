function [x, RESNORM, RESIDUAL, OUTPUT] = GaussNewton(rFun, x0,options)

% function [x, RESNORM, RESIDUAL, EXITFLAG, OUTPUT] = GaussNewton(rFun, x0, ...
%    options)
%
%	(c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function finds Gauss-Newton solution to an unconstrained optimization
% problem. The objective is to minimize the 2-norm of an m x 1 residual
% vector for an n x 1 optimization variable. The implementation follows
% notations on Page 209 of Lieven Vandenberghe's UCLA EE103 course reader.
%   http://www.seas.ucla.edu/~vandenbe/103/reader.pdf
%
% Input:
%   rFun: m x 1 objective residual vector r. The vector contains the m number
%       of residuals whose sum of squares are to be minimized.
%   x0: n x 1 initial starting point of the optimization variables
%   Jacobian: The m x n Jacobian matrix of the residuals with respect to
%       the optimzation variables (dr_i / dx_j)
%	options: user defined options for Gauss-Newton optimization process
%		maxNumItr: Maximum number of iterations for the Gauss-Newton
%		optimization; (default: 1e4)
%   	TolGrad: Termination tolerance on the gradient at current values
%				 of optimization variables; (default: 1e-3)
%   	TolX: Termination tolerance on the value change of the optimization
%      		  variables between two iterations; (default: 1e-5)
%   	TolFun: Termination tolerance on the value change of the objective
%      			function between two iterations;(default: 1e-5)
%       display:display intermediate optimization results;(default:on)
%		Constraint: lower and upper bound of the updating variables;(default: Inf)
%
% Output:
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
TolGrad = 1e-3;
TolX = 1e-5;
TolFun = 1e-5;
maxNumItr = 1e4;
lb = -Inf * ones(length(x0),1);
ub = Inf*ones(length(x0),1);
display = 'on';


%% Use user-defined options if exists
if(nargin == 3)
    if(isfield(options,'TolGrad'))
        TolGrad = options.TolGrad;
    end
    if(isfield(options,'TolX'))
        TolX = options.TolX;
    end
    
    if(isfield(options,'TolFun'))
        TolFun = options.TolFun;
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
    if(isfield(options,'display'))
        display = options.display;
    end
    
end

header = sprintf(['\n                                  Norm of      First-order \n',...
        ' Iteration          f(x)          step          optimality\n']);
if(strcmp(display,'on'))
    fprintf(header);
end

formatstr = ' %5.0f     %13.6g  %13.6g   %12.3g      \n';

alpha = 0.01; % slope factor

x = x0;

for k = 1 : maxNumItr
    [r,A] = feval(rFun,x);
    
    % g = sum of squares of the residuals r_i
    % Gradient of g over x: dg / dx
    grad = 2 * A' * r;
    if (norm(grad) < TolGrad)
        new_r = r;
        OUTPUT.exitflag = 1;
        break;
    end
    
    v = -A \ r;   % Decide the Gauss-Newton search direction
    t = ones(length(x),1);
        
    % Make nexX within the bounds
    temp1 =(x + t .* v < lb);
    temp2 =(x + t .* v > ub);
    while (~isempty(find(temp1,1)) || ~isempty(find(temp2,1)))
        lb_idx = find(temp1,1);
        ub_idx = find(temp2,1);
        t([lb_idx;ub_idx]) = t([lb_idx;ub_idx]) / 2;
        temp1 =(x + t .* v < lb);
        temp2 =(x + t .* v > ub);
        
    end
    new_x = x + t .* v;
    new_r = feval(rFun, new_x);
    
    % Backtracking step size
    while ( norm(new_r)^2 > norm(r)^2 + alpha * (grad' * v) * t )
        t = t/2;
        new_x = x + t .* v;
        new_r = feval(rFun, new_x);
    end
    x = x + t .* v;
    
    if(norm(t .* v) < TolX)
        OUTPUT.exitflag = 2;
        break;
    end
    
    if(abs(norm(new_r)^2 - norm(r)^2) < TolFun)
        OUTPUT.exitflag = 3;
        break;
    end
    if(strcmp(display,'on'))
        fprintf(formatstr,k,norm(new_r)^2,norm(t .* v),norm(grad))
    end
end

if (k == maxNumItr)
    OUTPUT.exitflag = 4;
end

RESNORM = norm(new_r)^2;
RESIDUAL = new_r;
OUTPUT.grad = grad;
OUTPUT.numIter = k;

end
