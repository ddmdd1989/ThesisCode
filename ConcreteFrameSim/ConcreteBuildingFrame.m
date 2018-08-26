clc;clear
close all
warning('off')


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = [-0.1;  0.2;   0.2; -0.05; 0.2;  0.15;
             0.15; 0.10; -0.10; -0.15; 0.20; 0.15;];

%% Load Sensitivity matrix
LoadStructure
n_modes = 3; % Number of measured modes
modeIndex = 1 : n_modes; % Indexes of these measured modes

%% Assemble structure matrices
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Optimization structure parameter;
optimzOpts.tolFun = 1e-10;
optimzOpts.tolX = 1e-10;
optimzOpts.toolBox = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e4;
optimzOpts.maxFunEvals = 2e4;


%% Simulate "experimental data"
[psiExp, lambdaExp] = eigs(K_act, M0, n_modes, 'sm');

[lambdaExp,dummyInd] = sort(diag(lambdaExp), 'ascend') ;
lambdaExp = lambdaExp(modeIndex);
psiExp = psiExp(:, dummyInd(modeIndex));
psi_m = psiExp(measDOFs,:);
unmeasDOFs = setdiff(1 : N, measDOFs);
psi_u  = psiExp(unmeasDOFs,:);

num_unmeasDOFs = N - length(measDOFs);
% Normalize the mode shape vectors by maximum entry
for i = 1: n_modes
    [~,index] = max(abs(psi_m(:,i)));
    psi_u(:,i) = psi_u(:,i) / psi_m(index,i);
    psi_m(:,i) = psi_m(:,i) / psi_m(index,i);
end

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
expModes.lambdaWeights = ones(n_modes,1);
expModes.psiWeights = ones(n_modes,1);
expModes.resWeights = ones(n_modes,1);

%% Model updating parameter
updatingOpts.formID = 3.0;       % 1: Modal property diff (MAC) ;
                                 % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID > 2.0 && updatingOpts.formID < 4.0)
    % Optimizaiton variable for modal dynamic residual formulation
    updatingOpts.x_lb = [-ones(n_alpha,1);-2 * ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub = [ones(n_alpha,1);  2 * ones(num_unmeasDOFs * n_modes,1)];
    
else
    % Optimizaiton variable for modal property difference formulations
    updatingOpts.x_lb = -ones(n_alpha,1);
    updatingOpts.x_ub =  ones(n_alpha,1);
end


%% MultiStart optimization
numRuns = 100;
randSeed = 3;
if(strcmp(optimzOpts.optAlgorithm,'Levenberg-Marquardt'))
    filename = ['ConcBuildFrm_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_LM.mat'];
else
    filename = ['ConcBuildFrm_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_TRR.mat'];
end

MultiRunModelUpdating


