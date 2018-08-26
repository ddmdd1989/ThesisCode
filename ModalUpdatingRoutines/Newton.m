function [x, RESNORM, RESIDUAL, OUTPUT] = Newton(rFun, x0, Jacobian, ...
    Hessian, options)

% function [x, RESNORM, RESIDUAL, EXITFLAG, OUTPUT] = GaussNewton(rFun, x0, ...
%    Jacobian, options)
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function finds Gauss-Newton solution to an unconstrained optimization
% problem. The objective is to minimize the 2-norm of an m x 1 residual
% vector for an n x 1 optimization variable. The implementation follows
% notations on Page 201 of Lieven Vandenberghe's UCLA EE103 course reader.
%   http://www.seas.ucla.edu/~vandenbe/103/reader.pdf
%
% Input:
%   rFun: function handle for scalar objective r.
%   x0: n x 1 initial starting point of the optimization variables
%   Jacobian: The n x 1 Jacobian vector of the objective function value
%             with respect to the optimization variables
%   Hessian: The n x n Hessian matrix of the objective function value
%             with respect to the optimization variables
%   options: user defined options for Newton optimization process
%       maxNumItr: Maximum number of iterations for the Gauss-Newton
%       optimization; (default: 1e4)
%       TolGrad: Termination tolerance on the gradient at current values
%                of optimization variables; (default: 1e-3)
%       TolX: Termination tolerance on the value change of the optimization
%             variables between two iterations; (default: 1e-5)
%       TolFun: Termination tolerance on the value change of the objective
%               function between two iterations;(default: 1e-5)
%
% Output:
%   x: optimal values of optimization variables
%   RESNORM: 2-norm of the optimal residual vector
%   RESIDUAL: optimal residual vector
%   OUTPUT: additional information of the optimization solution
%       grad: gradient of the sum of squares of r_i over x
%       numIter: Number of iterations of Gauss-Newton optimization process.
%       exitflag: exit conditions with following four values
%           1: TolGrad termination criterion
%           2: TolX termination criterion
%           3: TolFun termination criterion
%           4: Exceed the specified maximum number of iterations (maxNumItr)

%% Options default values
TolGrad = 1e-3;
TolX = 1e-5;
TolFun = 1e-5;
maxNumItr = 1e4;


%% Use user-defined options if exists
if(nargin == 4)
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
end


alpha = 0.01; % slope factor

x = x0;

for k = 1 : maxNumItr
    r = feval(rFun,x);
    J = feval(Jacobian,x);
    H = feval(Hessian,x);
    if(size(J,1) == 1)
        J = J';
    end
    
    if (norm(J) < TolGrad)
        new_r = r;
        OUTPUT.exitflag = 1;
        break;
    end;
    
    [~,p] = chol(H);
    if(p == 0)
        v = -inv(H) * J;
    else
        v = -J;
    end
    t = 1;
    new_x = x + t*v;
    new_r = feval(rFun, new_x);
    % Backtracking step size
    while ( new_r > r + alpha * t * J' * v)
        
        t = t/2;
        new_x = x + t * v;
        new_r = feval(rFun, new_x);
    end;
    x = x + t*v;
    
    if(norm(t*v) < TolX)
        OUTPUT.exitflag = 2;
        break;
    end
    
    if(abs(norm(new_r)^2 - norm(r)^2) < TolFun)
        OUTPUT.exitflag = 3;
        break;
    end
end

if (k == maxNumItr)
    OUTPUT.exitflag = 4;
end

RESNORM = norm(new_r);
RESIDUAL = new_r;
OUTPUT.grad = J;
OUTPUT.numIter = k;

end
