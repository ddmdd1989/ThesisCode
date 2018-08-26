function [J] = JacMDR(x,structModel,expModes)

% function [J] = JacMDR(x,M,K0,Lambda_e,V_m,P_K,MDOF,UDOF)
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function calculates the Jacobian matrix of the objective residual vector r
% for the objective function of the modal dynamic residual approach with respect to
% the optimization variables
% Note: the optimization variables include the stiffness parameters and the unmeasured
% entries of the eigenvectors
%
% Input:
%   x: the values of the optimization variables where the Jacobian matrix is calculated
%   M: mass matrix
%   K0: initial stiffness matrix
%   Lambda_e: experimental measured eigenvalues
%   V_m: experimental measured eigenvectors at measured DOF
%   P_K: influence matrices of the stiffness parameters
%   MDOF: measured DOFs
%   UDOF: unmeasured DOFs
% Output:
%   J: Jacobian matrix of the objective function

weight = (expModes.lambdaWeights + expModes.psiWeights) / 2;
N = size(structModel.M0,1);
n_alpha = size(structModel.K_j,3);
alpha = x(1:n_alpha);
n_u = N - expModes.n_meas; 
V_ep = zeros(N,expModes.n_modes);
V_u = zeros(n_u,expModes.n_modes);
for ii = 1 : expModes.n_modes
    V_u(:,ii) = x(n_alpha + (ii-1) * n_u + 1: n_alpha + ii * n_u);
    V_ep(1:expModes.n_meas,ii) = expModes.psiExp(:,ii);
    V_ep(expModes.n_meas + 1:end,ii) = V_u(:,ii);
end
for ii = 1 : expModes.n_modes
    S_u0 = zeros(N,n_u);
    for jj = 1:n_alpha
        S_u0 = S_u0 + alpha(jj) * structModel.K_j(:,expModes.n_meas + 1:end,jj);
    end
    J_v((ii-1) * N + 1 : ii * N , (ii-1) * n_u + 1: ii * n_u) = (structModel.K0(:,expModes.n_meas + 1:end) + S_u0...
        - expModes.lambdaExp(ii) * structModel.M0(:,expModes.n_meas + 1:end)) * weight(ii);
end
for ii = 1:n_alpha
    for jj = 1:expModes.n_modes
        J_alpha((jj-1) * N + 1:jj * N,ii) = structModel.K_j(:,:,ii) * V_ep(:,jj) * weight(jj);
    end
end

J = [J_alpha J_v];