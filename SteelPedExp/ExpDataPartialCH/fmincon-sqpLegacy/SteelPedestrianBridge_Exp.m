%% 
clear
close all

LoadStructure;

n_alpha = size(K_j,3) - 2;
modeIndex = [1 2 4];

%% Assemble structure parameter:
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j(:,:,1 : n_alpha);

%% Optimization structure parameter;
optimzOpts.tolFun = 1e-10;
optimzOpts.tolX = 1e-10;
optimzOpts.tolGrad = 1e-10;
optimzOpts.toolBox = 'fmincon';
% optimzOpts.optAlgorithm = 'trust-region-reflective';
optimzOpts.optAlgorithm = 'sqp-legacy';

optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 1e3;
optimzOpts.maxFunEvals = 3e3;

%% Model updating parameter
updatingOpts.formID = 1;       % 1: Modal property diff (MAC) ;
                               % 2: Modal property diff (V_mDiff);
updatingOpts.modeMatch = 2;   % 1: Without forced matching;
                               % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;

updatingOpts.x_lb = [-1* ones(n_alpha - 3,1); -1;   -ones(2,1)];
updatingOpts.x_ub = [ 1 * ones(n_alpha - 3,1);  1;  10 * ones(2,1) ];

updatingOpts.reg = 15;


%% Experimental modal properties
load ModeInfo

lambdaExp = Lambda_mean(modeIndex);

psiExp_m = mode_mean(:,modeIndex);

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psiExp_m;
expModes.measDOFs = measDOFs;
expModes.good_ch = good_ch(modeIndex);

for i = 1:length(modeIndex)
    if(updatingOpts.formID == 4 || updatingOpts.formID == 1)
        weightPsi(i) = 1 / norm(mode_std(good_ch{modeIndex(i)},modeIndex(i)));
    else
        weightPsi(:,i) = 1 ./ mode_std(:,modeIndex(i));
    end
end
weightLambda = lambdaExp ./ std_Lambda(modeIndex);
expModes.lambdaWeights = weightLambda .* [2;1;3];
expModes.psiWeights = weightPsi;

if(strcmp(optimzOpts.toolBox,'L_M'))
    filename = ['Space_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_LMOwn.mat'];
elseif(strcmp(optimzOpts.toolBox,'Gauss-Newton'))
    filename = ['Space_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_GSNT.mat'];
elseif(strcmp(optimzOpts.toolBox,'lsqnonlin'))
    if(strcmp(optimzOpts.optAlgorithm,'Levenberg-Marquardt'))
        filename = ['Space_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_LM.mat'];
    elseif(strcmp(optimzOpts.optAlgorithm,'trust-region-reflective'))
        filename = ['Space_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_TRR.mat'];
    end
elseif(strcmp(optimzOpts.toolBox,'fmincon'))
    filename = ['Space_form' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_intPoint.mat'];
end

numRuns = 100 ;
randSeed = 1;
MultiRunModelUpdating
% 

modeIndex = [1 2 4];
unmeasDOFs = setdiff(1 : N, measDOFs);
rordIdx = [measDOFs';unmeasDOFs'];
K0 = K0(rordIdx, rordIdx);
M0 = M0(rordIdx,rordIdx);
K_j = K_j(rordIdx,rordIdx,:);
f_e = sqrt(lambdaExp) / 2 / pi;


[psiNom,lambdaNom] = eigs(K0, M0, 20, 'sm');
[lambdaNom, srt_idx] = sort(diag(lambdaNom),'ascend');
psiNom_m = psiNom(1 : length(measDOFs), srt_idx);
psiNom_m = psiNom_m(:, modeIndex);
lambdaNom = lambdaNom(modeIndex);

f_nom = sqrt(lambdaNom) / 2/ pi;

MAC_nom(1) = mac(psiExp_m(good_ch{1},1),psiNom_m(good_ch{1},1));
MAC_nom(2) = mac(psiExp_m(good_ch{2},2),psiNom_m(good_ch{2},2));
MAC_nom(3) = mac(psiExp_m(good_ch{4},3),psiNom_m(good_ch{4},3));

[~,idx] = min(fval);
opt = alpha(:,idx);

K_opt = sum(K_j,3);

for i = 1:n_alpha
    K_opt = K_opt + K_j(:,:,i) * opt(i);
end

[psiOpt,lambdaOpt] = eigs(K_opt,M0,20,'sm');
[lambdaOpt,srt_idx] = sort(diag(lambdaOpt),'ascend');
psiOpt_m = psiOpt(1 : length(measDOFs), srt_idx);
psiOpt_m = psiOpt_m(:,modeIndex);
lambdaOpt = lambdaOpt(modeIndex);
f_opt = sqrt(lambdaOpt) / 2/ pi;

MAC_opt(1) = mac(psiExp_m(good_ch{1},1),psiOpt_m(good_ch{1},1));
MAC_opt(2) = mac(psiExp_m(good_ch{2},2),psiOpt_m(good_ch{2},2));
MAC_opt(3) = mac(psiExp_m(good_ch{4},3),psiOpt_m(good_ch{4},3));


