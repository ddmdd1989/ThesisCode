clear
close all
clc



%% Strucutural parameters
n_modes = 3;  % number of measured mode shapes
N = 6;   % DOF of the whole structure
masses = 6 * ones(N, 1);             % kg
iniSpring = 35000 * ones(N, 1);          % N/m

%% Simulated actual structure
dmgLoc = 1:6;

alpha_act = [-0.25; -0.20; -0.15; -0.10; -0.05; 0.05];

actSpring = iniSpring;

for i = 1:length(dmgLoc)
    actSpring(dmgLoc(i)) = actSpring(dmgLoc(i)) * (1 + alpha_act(i));
end


%% Actual strucuture
M0 = makeM(masses, N);
K_act = makeK(actSpring, N);
[psiExp,lambdaExp] = eigs(K_act, M0,5,'sm') ;
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;
measDOFs = [1 3 6];
unmeasDOFs = setdiff(1:N,measDOFs);
num_unmeasDOFs = length(unmeasDOFs);
lambdaExp = lambdaExp(1:n_modes);
modeIndex = 1 : n_modes;
psiExp = psiExp(:,dummyInd);

psiExp_m = psiExp(measDOFs,1:n_modes) ;     % simulated measurement
psiExp_u = psiExp(unmeasDOFs,1:n_modes) ;    
for i = 1:n_modes
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    psiExp_u(:,i) = psiExp_u(:,i) / psiExp_m(q(i),i);
    psiExp_m(:,i) = psiExp_m(:,i) / psiExp_m(q(i),i);
end

%% Initial strucuture modes
M0 = makeM(masses, N);
K0 = makeK(iniSpring, N);

%% Assemble sensitivty matrix
i = 1 ;
K_j(:,:,i) = zeros(N) ;
K_j(i,i,i) = iniSpring(i);

for i = 2 : N
    K_j(:,:,i) = zeros(N) ;
    K_j(i-1,i-1,i) = iniSpring(i);
    K_j(i-1,i,i) = -iniSpring(i);
    K_j(i,i-1,i) = -iniSpring(i);
    K_j(i,i,i) = iniSpring(i);
end

n_alpha = size(K_j,3);


%% Assemble structure matrices
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j;

%% Optimization structure parameter;
optimzOpts.tolFun = eps^2;
optimzOpts.tolX = eps^2;
optimzOpts.toolBox = 'lsqnonlin';
optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e3;
optimzOpts.maxFunEvals = 3e5;


%% Simulate "experimental data"
expModes.lambdaExp = lambdaExp;
expModes.psiExp = psiExp_m;
expModes.measDOFs = measDOFs;
expModes.lambdaWeights = ones(n_modes,1);
expModes.psiWeights = ones(n_modes,1);
expModes.resWeights = ones(n_modes,1);

%% Model updating parameter
updatingOpts.formID = 3.0;       % 1: Modal property diff (MAC) ;
                                 % 2: Modal property diff (V_mDiff);
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
randSeed = 2;
if(strcmp(optimzOpts.optAlgorithm,'Levenberg-Marquardt'))
    filename = ['LumpedSpringMass_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_LM.mat'];
else
    filename = ['LumpedSpringMass_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_TRR.mat'];
end

MultiRunModelUpdating




