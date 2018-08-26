function [c,ceq]= FminconCons(x,K0,M,P_K,MDOF,UDOF,Lambda_e,V_m,objOpt)

% function [c,ceq]= ObjFuncMPDFminconCons(x,K0,M,P_K,MDOF,UDOF,Lambda_e,V_m,,objOpt)
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function calculates the nonlinear equality constraints included in the optimization problems
% of the modal property difference approach. 
% Note:
%   1. The function only intends to update the stiffness parameters of the structure
%   2. The optimization variables included:
%       a. structure parameters
%       b. analytical eigenvalues 
%       c. measured entries of the analytical eigenvectors
%       d. unmeasured entries of the analytical eigenvectors
%       e. infinity norm of the modal property difference vector
%   3. The function is used when the optimization problem of the modal property difference approach 
%      be applied with fmincon in MATLAB 
% Input:
%   x: optimization variables
%   K0: initial stiffness matrix
%   M: mass matrix
%   P_K: influence matrices of updating parameters
%   MDOF: measured DOF of eigenvectors
%   UDOF: unmeasured DOF of eigenvectors
%   Lambda_e: experimental eigenvalues
%   V_m: experimental eigenvector at measured DOF
%   objOpt:
%   1. feasibility problem, set modal property difference vector = 0;
%   2. minimize the infinity norm of the modal property difference vector;
% Output:
%   c: the nonlinear inequality constraints involved in the optimization problem of the modal property
%        difference approach
%   ceq: the nonlinear equality constraints involved in the optimization problem of the modal property
%        difference approach

n_alpha = size(P_K,3);
numModes = length(Lambda_e) ;
nMDOF = length(MDOF);
nUDOF = length(UDOF);
N = nMDOF + nUDOF;
K = K0;
V_a = zeros(N,1);
n_m = size(V_m,1);
for i = 1:n_alpha
    K = K + x(i) * P_K(:,:,i);
end
for i = 1:numModes
    Lambda_a = x(n_alpha + i);
    V_a(MDOF,1) = x(n_alpha + numModes + (i-1) * nMDOF + 1:n_alpha + numModes + i * nMDOF);
    V_a(UDOF,1) = x(n_alpha + numModes + numModes *  nMDOF + (i-1) * nUDOF  + 1 : n_alpha + numModes + numModes *  nMDOF + i * nUDOF);
    y((i-1) *N + 1: i * N,1)  = ( K - Lambda_a * M) * V_a;
    
end
for i = 1:numModes
    Lambda_a(i,1) = x(n_alpha + i);
    V_ma(:,i) = reshape(x(n_alpha + numModes + (i-1) * n_m + 1 : n_alpha + numModes + i * n_m),n_m,1);
end

MACValue = diag(mac(V_ma,V_m)) ;
obj = [(Lambda_a - Lambda_e)./Lambda_e;  (1-sqrt(MACValue))./ sqrt(MACValue)];

if(objOpt == 1)
    c =0;
    ceq = [y;obj];
else
    index = n_alpha + numModes + numModes * N;
    c = [obj - x(index + 1) * ones(length(obj),1);-obj - x(index + 1) * ones(length(obj),1)];
    ceq = y;
    
end
end




