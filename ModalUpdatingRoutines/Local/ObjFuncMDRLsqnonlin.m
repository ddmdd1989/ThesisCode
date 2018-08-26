function r = ObjFuncMDRLsqnonlin(x,structModel,expModes)

% function [r] = ObjFuncMDRLsqnonlin(alpha,structModel,expModes,opt)
%
%   (c) Yang Wang, Dapeng Zhu, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function calculates the objective function of the modal dynamic 
%   residual approach 
% Note: this function only intends to update the stiffness matrix
% 
% Input:
%   x: optimization variables, including stiffness parameters and 
%       unmeasured entries of the eigenvectors
%   alpha - optimization variables
%   structModel - a structure array with following fields of structural
%   model information:
%       M0 (N x N)- mass matrix (assumed accurate enough and no need to
%          update in current revision). Here N refers to the number of
%          degrees of freedom of the finite element model
%       K0 (N x N) - nominal stiffness matrix constructed with nominal
%          parameter values
%       K_j (N x N x n_alpha) - influence matrix corresponding to updating
%          variables (Note: the third dimension of K_j should be
%          equal to the number of updating variables). Here n_alpha refers
%          the number of stiffness updating variables
%
%   expModes - a structure array with experimental modal properties for
%   model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs.Here n_meas refers to the number of measured DOFs 
%       measDOFs (n_meas x 1) - measured DOFs  
%       lambdaWeights (n_modes x 1) - weighting factor for eigenvalue
%       psiWeights (n_modes x 1) - weighting factor for eigenvector
% Output:
%   y: the output of the objective function(vector)

weight = expModes.resWeights;
N = size(structModel.K_j,1);
n_modes = length(expModes.lambdaExp);
psiMix = zeros(N,1);

n_alpha = size(structModel.K_j,3);

n_u = N - expModes.n_meas; 
r = zeros( N * n_modes,1);
for i = 1 : n_modes
    psiMix(1:expModes.n_meas,1) = expModes.psiExp(:,i);
    psiMix(expModes.n_meas + 1 : end,1) = x(n_alpha + (i-1) * n_u + 1: n_alpha + i * n_u);
    r((i-1) * N + 1 : i * N,:) = (structModel.K - structModel.M0 * expModes.lambdaExp(i))* psiMix * weight(i);
end

end

