function r = ModelUpdatingObjective(x, structModel, expModes, ...
    simModes, updatingOpts)
% function r = ModelUpdatingObjective(alpha, structModel, expModes,
%   simModes, updatingOpts)
%
%   (c) Yang Wang, Xinjun Dong, Dan Li (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% This function calculates the objective residuals of the
% optimization problem for model updating. When MATLAB lsqnonlin solver is
% used, the output contains a vector (r) whose entries are the residuals.
%
% Input:
%   x - a vector with values of the optimization variables
%
%   structModel - a MATLAB structure array with following fields of
%   structural model information:
%       M0 (N x N)- mass matrix (assumed accurate enough and no need to
%          update in current revision). Here N refers to the number of
%          degrees of freedom of the finite element model
%       K0 (N x N) - nominal stiffness matrix constructed with nominal
%          parameter values
%       K_j (N x N x n_alpha) - influence matrix corresponding to updating
%          variables (Note: the third dimension of K_j should be
%          equal to the number of updating variables). Here n_alpha refers
%          the number of stiffness updating variables
%       K (N x N) - stiffness matrix constructed with the current alpha
%         values, using K0 and K_j
%
%   expModes - a MATLAB structure array with experimental modal properties
%   for model updating:
%       lambdaExp (n_modes x 1) - experimental eigenvalue. Here n_modes
%          refers to the number of experimental modes available
%       psiExp (n_meas x n_modes) - experimental mode shape vector at
%          measured DOFs. Here n_meas refers to the number of measured DOFs
%       measDOFs (n_meas x 1) - measured DOFs
%       lambdaWeights (n_modes x 1) - weighting factor for eigenvalue
%       psiWeights (n_modes x 1) - weighting factor for eigenvector
%
%   simModes - a structure array with simulated modal properties for
%     model updating:
%       Lambda (n_modes x 1) - simulated eigenvalue
%       psi_m  (n_meas x n_modes) - simulated mode shape vector at
%          measured DOFs
%       psi    (N x n_modes) - simulated mode shape vector at all DOFs
%
%   updatingOpts - a MATLAB structure array with model updating options:
%       formID - formulation ID number (default: 1)
%           1: Case 1 - conventional modal property difference
%              formulation using MAC values
%           2: Case 2 - modal property difference formulation with
%             eigenvector difference formulation
%
% Output:
%	r: the objective residual vector r(x)

switch updatingOpts.formID
    case 1.0
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 0, 1, 1);
    case 1.1
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 1, 1, 1);
    case 1.2
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 2, 1, 1);
    case 2.0
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 0, 1, 2);
    case 2.1
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 1, 1, 2);
    case 2.2
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 2, 1, 2);
    case 2.3
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 0, 2, 2);
    case 2.4
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 1, 2, 2);
    case 2.5
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 2, 2, 2);
    case 3.0
        r =  ObjFuncMDRLsqnonlin(x, structModel, expModes, 1);
    case 4.0
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 0, 1, 1);
    case 4.1
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 1, 1, 1);
    case 4.2
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 2, 1, 1);
    case 5.0
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 0, 1, 2);
    case 5.1
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 1, 1, 2);
    case 5.2
        r = ObjFuncMPDLsqnonlin(expModes, simModes, 2, 1, 2);
    otherwise
        error('\nWrong option for objective function!');
end

end
