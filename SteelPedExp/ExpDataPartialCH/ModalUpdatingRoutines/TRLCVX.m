function [x, RESNORM, OUTPUT] = TRLCVX(rFun,Jacobian,Hessian,x0,options)

% function [x, RESNORM, RESIDUAL, EXITFLAG, OUTPUT] = TRLCVX(rFun,Jacobian,Hessian,x0,options)
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function finds Trust-Region-Reflective solution to an unconstrained optimization
% problem. The objective is to minimize the 2-norm of an m x 1 residual
% vector for an n x 1 optimization variable.
%
% Input:
%   rFun: an objective function returns a scalar value which is to be minimized.
%   x0: n x 1 initial starting point of the optimization variables
%   Jacobian: The n x 1 Jacobian vector of the scalar with respect to
%       the optimization variables (dr / dx_j)
%   options: user defined options for Gauss-Newton optimization process
%       maxNumItr: Maximum number of iterations for the Gauss-Newton
%       optimization; (default: 1e4)
%       TolGrad: Termination tolerance on the gradient at current values
%                of optimization variables; (default: 1e-3)
%       TolX: Termination tolerance on the value change of the optimization
%             variables between two iterations; (default: 1e-5)
%
% Output:
%   x: optimal values of optimization variables
%   RESNORM: the optimal objective function value
%   OUTPUT: additional information of the optimization solution
%       grad: the optimal gradient value
%       numIter: Number of iterations of Trust-Region-Reflective optimization process.
%       exitflag: exit conditions with following four values
%           1: TolGrad termination criterion
%           2: TolX termination criterion           
%           3: TolFun termination criterion
%           4: Exceed the specified maximum number of iterations (maxNumItr)


%% Options default values
TolGrad = 1e-3;
TolX = 1e-5;

maxNumItr = 1e4;

%% Use user-defined options if exists
if(nargin == 5)
    if(isfield(options,'TolGrad'))
        TolGrad = options.TolGrad;
    end
    if(isfield(options,'TolX'))
        TolX = options.TolX;
    end
    
    if(isfield(options,'maxNumItr'))
        maxNumItr = options.maxNumItr;
    end
end




RRT  = 1/5;        % Reduction Ratio Threshhold ("eta" in notes).
RRT_L = 0.25;      % Lower bound of reduction ratio, otherwise shrink trust region
RRT_U = 0.75;      % Upper bound of reduction ratio, otherwise increase trust region
Delta = 0.01;       % Initial value for trust region radius.
Dmax = 1.5;        % Absolute maximum trust region radius.


n_alpha = length(x0);

xk = x0;

for k = 1:maxNumItr,
    %% Calculate the stepsize within the assigned bound
    % Cauchy point method
    % Evaluate the F,J and H
    
    fk = feval(rFun,xk);
    gk = feval(Jacobian,xk);
    Hk = feval(Hessian,xk);
    
    
    [U,D]=eig(Hk);
    f=U'*gk';
    d=diag(D);
    
    cvx_begin quiet
    variable z(n_alpha)
    minimize(d'*z-2*abs(f)'*sqrt(z))
    subject to
    sum(z) <= Delta
    z >= 0
    cvx_end
    y=-sign(f).*sqrt(z);
    P_c=U*y;
    
    
    
    % Estimated function reduction
    predred = -(gk * P_c + 0.5 * P_c' * Hk * P_c);
    % Actual function reduction
    xc = xk + P_c;
    ft = feval(rFun,xc);
    
    actred = fk - ft;   % Actual function reduction
    
    RR = actred/predred;   % Reduction Ratio.
    
    
    % Update the current point
    if RR >= RRT,
        xk = xc;
    end
    
    % Update the trust region radius:
    if RR < RRT_L
        Delta = 0.25 * Delta;
    elseif RR > RRT_U
        Delta = min([2*Delta;Dmax]);
    end
    
    
    
    if (norm(P_c) < TolX)
        OUTPUT.exitflag = 2;
        break;
    end
    if (norm(gk) < TolGrad)
        OUTPUT.exitflag = 1;
        break;
    end
   
    if (k == maxNumItr)
        OUTPUT.exitflag = 4;
    end
    
    display(ft);
end

x = xk;
RESNORM = ft;
OUTPUT.grad = gk;
OUTPUT.numIter = k;
end




