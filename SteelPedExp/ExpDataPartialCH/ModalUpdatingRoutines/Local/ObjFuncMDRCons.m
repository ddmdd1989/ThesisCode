function [c,ceq] = ObjFuncMDRCons(x,K0,M,P,Lambda_e,V_m,MDOF,UDOF,weight,objOpt)

% function [c,ceq]= ObjFuncMDRCons(x,K0,M,P,Lambda_e,V_m,MDOF,UDOF,weight,objOpt)
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function calculates the nonlinear equality constraints included in 
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
%   K0: initial stiffness matrix
%   M: mass matrix
%   P: influence matrices of updating parameters
%   Lambda_e: experimental eigenvalues
%   V_m: experimental eigenvector at measured DOF
%   MDOF: measured DOF of eigenvectors
%   UDOF: unmeasured DOF of eigenvectors
%   weight: weighing factors for the measured modes
%   objOpt:
%   1. feasibility problem, set modal dynamic residual vector = 0;
%   2. minimize the infinity norm of the modal dynamic residual vector;
% Output:
%   c: the nonlinear inequality constraints between the modal dynamic residual and the upper bounds of  
%      the infinite norm of the modal dynamic residual vector
%   ceq: the nonlinear equality constraints involved in the optimization problem of the modal dynamic
%        residual approach


N_tilde = size(M,1) ;
numModes = size(V_m,2) ;
V_ep = zeros(N_tilde,numModes);
K = K0;
n_alpha = size(P,3);
n_u = length(UDOF);

for i = 1 : n_alpha
    K = K + x(i) * P(:,:,i) ;
end

for i = 1 : numModes
    V_ep(MDOF,i) = V_m(:,i) ;
    V_ep(UDOF,i) = x(n_alpha + (i-1) * n_u+1 : n_alpha + i * n_u)';
end
% index = n_alpha + numModes * length(UDOF);
for j = 1 : numModes
    y((j-1)*N_tilde+1:j*N_tilde,1) = (K-Lambda_e(j)*M)*V_ep(:,j) * weight(j);
end

if(objOpt == 1)
    c = 0;
    ceq = y;
else
    index = n_alpha + numModes * length(UDOF);
    c = [y - x(index + 1) * ones(length(y),1);-y - x(index + 1) * ones(length(y),1)];
    ceq = 0;
end
end
