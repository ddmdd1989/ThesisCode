function [x, RESNORM, OUTPUT] = TRLCVXBoxCons(rFun,Jacobian,Hessian,x0,options)

% function [x, RESNORM, RESIDUAL, EXITFLAG, OUTPUT] =
%   TRLCVXBoxCons(rFun,Jacobian,Hessian,x0,options) 
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function finds Trust-Region-Reflective solution to an optimization
% problem with box constraints. The objective is to minimize the 2-norm of
% an m x 1 residual vector for an n x 1 optimization variable.
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
%       display:display intermediate optimization results;(default:on)
%       Constraint: lower and upper bound of the updating variables;(default: Inf)
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
display = 'on';
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
    
    xc = xk + P_c;
    
        % Make nexX within the bounds
    temp1 =(xk + P_c < lb);
    temp2 =(xk + P_c > ub);
    while (~isempty(find(temp1,1)) || ~isempty(find(temp2,1)))
        P_c = P_c / 2;
        temp1 =(xk + P_c < lb);
        temp2 =(xk + P_c > ub);
        
    end;
    xc = xk + P_c;
      
    
    fc = feval(rFun,xc);
    
    actred = fk - fc;   % Actual function reduction
    
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
    if(strcmp(display,'on'))
        fprintf(formatstr,k,fc,norm(P_c),norm(gk))
    end
    
end

x = xk;
RESNORM = fc;
OUTPUT.grad = gk;
OUTPUT.numIter = k;
end




