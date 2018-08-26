function y = ObjFuncMPDLsqnonlinModalMass(x,M0,K0,Lambda_e,V_m,P_K,P_M,MDOF,weight,opt)
% function y = ObjFuncMPDLsqnonlinModalMass(x,M0,K0,Lambda_e,V_m,P_K,P_M,MDOF,weight,opt)
%
%   (c) Yang Wang, Dapeng Zhu, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function calculates the objective function of the modal property difference approach
% Note: this function intends to update the stiffness and modal mass of the structure; modal
% mass is only used in substructure model updating
%
% Input:
%   x: optimization variables
%   M: mass matrix
%   K0: initial stiffness matrix
%   Lambda_e: experimental eigenvalues
%   V_m: experimental eigenvector at measured DOF
%   P_K: influence matrices of updating parameters
%   MDOF: measured DOF of eigenvectors
%   weight: weighing factors of the measured modes
%   opt: output format of the objective function
%       1. modal property difference vector
%       2. 2-norm of the modal property difference vector
% Output:
%   y: the output of the objective function of modal property difference approach (vector/scalar)
%      EWSHM 2016      Section 2.4

K = K0; M = M0;

n_alpha = size(P_K,3) ;
n_beta = size (P_M,3) ;
m = length(freq_e) ;

for i = 1 : n_alpha
    K = K + x(i) * P_K(:,:,i) ;
end

for i = 1 : n_beta
    M = M + x(i+n_alpha) * P_M(:,:,i) ;
end

[V_a,Lambda_a] = eig(K, M ) ;
[Lambda_a,I] = sort((diag(Lambda_a)),'ascend') ;
Lambda_a = Lambda_a(1:m);
V_a = V_a(:,I);

MACValue = diag(mac(V_a(MDOF,1:m),V_m)) ;

y = [((Lambda_a - Lambda_e)./Lambda_e) .*weight;   (1-sqrt(MACValue))./ sqrt(MACValue).*weight];
if(opt == 2)
    y = norm(y);
end

