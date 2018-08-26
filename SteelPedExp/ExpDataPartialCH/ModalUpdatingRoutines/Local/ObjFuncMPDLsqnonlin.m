function [r] =  ObjFuncMPDLsqnonlin(expModes, simModes, eigFreqOpt, normOpt, objOpt)
% function [r] =  ObjFuncMPDLsqnonlin(expModes, simModes, eigFreqOpt, normOpt, objOpt)
%
%   (c) Yang Wang, Xinjun Dong, Dan Li (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% For implementation with MATLAB lsqnonlin, this function calculates
% the objective residual vector for various forms of the modal property
% difference approaches.
%
% Input:
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
%   simModes - a MATLAB structure array with simulated modal properties for
%     model updating:
%       Lambda (n_modes x 1) - simulated eigenvalue
%       psi_m  (n_meas x n_modes) - simulated mode shape vector at
%          measured DOFs
%       psi    (N x n_modes) - simulated mode shape vector at all DOFs
%
%   eigFreqOpt
%       0 - use eigenvalue difference
%       1 - angular frequency difference (rad/s)
%       2 - ordinary frequency differnce (Hz)
%
%   normOpt
%       1 - normalize the qi-th entry of the eigenvector to 1
%       2 - normalize the length of the eigenvector to 1
%
%   objOpt
%       1 - MAC value formulation
%       2 - eigenvector difference formulation

% Output:
%	r: the objective residual vector r(x)

omegaSim = sqrt( simModes.Lambda );  % simulated angular frequency
omegaExp = sqrt( expModes.lambdaExp );% experimental angular frequency
freqSim = omegaSim / 2 / pi;
freqExp = omegaExp / 2 / pi;

n_meas = expModes.n_meas;
n_modes = expModes.n_modes;
r = [];
for i = 1 : n_modes
    if (eigFreqOpt == 0)
        eigFreqTerm = expModes.lambdaWeights(i) * ...
            (expModes.lambdaExp(i) - simModes.Lambda(i)) / expModes.lambdaExp(i);
    elseif (eigFreqOpt == 1)
        eigFreqTerm = expModes.lambdaWeights(i) * (omegaExp(i) - omegaSim(i)) / omegaExp(i);
    elseif (eigFreqOpt == 2)
        eigFreqTerm = expModes.lambdaWeights(i) * (freqExp(i) - freqSim(i)) / freqExp(i);
    end
    
    if (objOpt == 1) % MAC value formulation
        psiExp = expModes.psiExp(expModes.good_ch{i},i);
        psiSim = simModes.psi_m(expModes.good_ch{i},i);
        MACi = MAC(psiExp, psiSim);
        MACterm = expModes.psiWeights(i) * (1 - sqrt(MACi)) ./ sqrt(MACi);
        r = [r;eigFreqTerm; MACterm];
        
    else % eigenvector difference formulation
        if (normOpt == 1) % Maximum entry is normalized to 1
            % This is a simplified implementation with constructing the
            % matrix Q_i in the formulation. Effectively this finds the
            % Q_i * {Psi^EXP - Psi^m} * weight}in the formulation
            psiExp = expModes.psiExp(expModes.good_ch{i},i);
            psiSim = simModes.psi_m(expModes.good_ch{i},i);
            psiWeight = expModes.psiWeights(expModes.good_ch{i},i);
            Q_i = setdiff(1:length(expModes.good_ch{i}), find(expModes.good_ch{i} == expModes.q(i)));
            psiDiff = (psiExp(Q_i) - psiSim(Q_i)) .* psiWeight(Q_i);
            r = [r; eigFreqTerm; psiDiff];
        else % Vector length is normalized to 1
        	% Ref: Ref: Vanik, Michael W. etc. "Bayesian probabilistic approach to structural health monitoring." 
        	% Journal of Engineering Mechanics 126.7 (2000): 738-745.
            innPrdct = expModes.psiExp(:,i)' * simModes.psi_m(:,i);
            psiDiff = (expModes.psiExp(:,i) - innPrdct * simModes.psi_m(:,i)) ...
                * expModes.psiWeights(i);
            r((i - 1) * (n_meas + 1) + 1 : i * (n_meas + 1)) = [eigFreqTerm; psiDiff];
        end
    end
end
