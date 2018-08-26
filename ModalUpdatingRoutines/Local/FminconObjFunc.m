function y = ObjFuncMPDFmincon(x,numModes,N,n_alpha,objOpt)

% function y = ObjFuncMPDFmincon(x,numModes,N,n_alpha,objOpt)
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function calculates the objective function value included in 
% the optimization problems of the modal dynamic residual approach. 
% Note:
%   1. The function only intends to update the stiffness parameters of the structure
%   2. The optimization variables included:
%       a. structure parameters
%       b. unmeasured entries of the analytical eigenvectors
%       c. infinity norm of the modal dynamic residual vector
%   3. The function is used when applying modal dynamic residual approach
%   with fmincon
% Input:
%   x: optimization variables
%   numModes: number of measured modes
%   N: number of DOF of the structure
%   n_alpha: number of stiffness parameters to be updated
%   objOpt:
%   1. feasibility problem, set modal dynamic residual vector = 0;
%   2. minimize the infinity norm of the modal dynamic residual vector;
% Output:
%   y: objective function value
if(objOpt == 1)
    y = 0;
else
    index = n_alpha + numModes + numModes * N;
    y  = x(index + 1);
end

end

