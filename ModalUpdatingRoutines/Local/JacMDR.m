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

weight = expModes.resWeights;
N = size(structModel.M0,1);
n_alpha = size(structModel.K_j,3);
n_u = N - expModes.n_meas;

Dr_alpha = zeros(N * expModes.n_modes, n_alpha);

for i = 1:n_alpha
    for j = 1:expModes.n_modes
        psiMix(1 : expModes.n_meas,1) = expModes.psiExp(:,j);
        psiMix(expModes.n_meas + 1 : N,1) = x(n_alpha + (j-1) * n_u + 1: n_alpha + j * n_u);
        Dr_alpha((j-1) * N + 1 : j * N, i) = structModel.K_j(:,:,i) * psiMix * weight(j);
    end
end

Dr_psiU = zeros(N * expModes.n_modes, n_u * expModes.n_modes);
for i = 1 : expModes.n_modes
    K_u = structModel.K(:, expModes.n_meas + 1: end);
    M_u = structModel.M0(:, expModes.n_meas + 1: end);
    Dr_psiU((i-1) * N + 1 : i * N , (i-1) * n_u + 1: i * n_u) = (K_u - expModes.lambdaExp(i) * M_u) * weight(i);
end



J = [Dr_alpha Dr_psiU];