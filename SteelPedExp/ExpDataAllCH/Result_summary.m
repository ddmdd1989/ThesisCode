%% 
clear
close all
cd ..
LoadStructure;

n_alpha = size(P_K,3) - 4;
Modeindex = [1 2 4];

%% Assemble structure parameter:
structModel.M0 = Morig;
structModel.K0 = K;
structModel.K_j = P_K(:,:,1:n_alpha);


%% Experimental modal properties
load ModeInfo0724



numModes = length(Modeindex);
numJoint = 46;
N = numJoint * 6 - 2  ;
dof = load('DOF.txt');

% 
measNodes = [36 22 1 4 38 24 7 8 9 10 40 26 11 12 42 28 13 14 15 16 ...
             44 30 19 20 46 21 35 33 34 32];

numdir = 2;
for i = 1 : length(measNodes)
    for j = 1 : numdir
        measDOFsOrig_wz((i-1)*numdir+j) = dof(measNodes(i),2+j);
    end
end

MDOF = measDOFsOrig_wz(measDOFsOrig_wz~=0) ;

Lambda_e = Lambda_mean(Modeindex);

V_m = mode_mean(:,Modeindex);

cd UpdatingResults_wSpring
% 
load Space_form2_JACon_TRR

UDOF = setdiff(1:N,MDOF);
reord = [MDOF';UDOF'];
K = K(reord,reord);
Morig = Morig(reord,reord);
P_K = P_K(reord,reord,:);
[V,D] = eigs(K,Morig,20,'sm');
[Lambda_a,srt_idx] = sort(diag(D),'ascend');


V_ma = V(1:length(MDOF),srt_idx);
V_ma = V_ma(:,Modeindex);
Lambda_a = Lambda_a(Modeindex);
f_e = sqrt(Lambda_e) / 2 / pi;
f_a = sqrt(Lambda_a) / 2/ pi;

eigvalDiffIni = (Lambda_a - Lambda_e) ./ Lambda_e .* expModes.lambdaWeights;
                            
eigvecDiffIni = [];
for i = 1:length(Modeindex)
    V_eG = V_m(expModes.good_ch{i},i);
    idx = find(V_eG == 1);
    V_mG = V_ma(expModes.good_ch{i},i);
    V_mG = V_mG / V_mG( idx);
    eigvecDiffIni = [eigvecDiffIni;(V_eG - V_mG) .* expModes.psiWeights(expModes.good_ch{i}, i)];
end

eigvecDiffIni = eigvecDiffIni(~isnan(eigvecDiffIni));
[~,idx] = min(fval);
opt = alpha(:,idx);


n_alpha = size( structModel.K_j, 3 );
K_opt = K;
for i = 1 : n_alpha
    K_opt = K_opt + opt(i) * P_K(:, :, i);
end

[V,D] = eigs(K_opt,Morig,20,'sm');
[Lambda_opt,srt_idx] = sort(diag(D),'ascend');
V_opt = V(1:length(MDOF),srt_idx);
V_opt = V_opt(:,Modeindex);
Lambda_opt = Lambda_opt(Modeindex);

eigvalDiffOpt = (Lambda_opt - Lambda_e) ./ Lambda_e .* expModes.lambdaWeights;

eigvecDiffOpt = [];
for i = 1:length(Modeindex)
    V_eG = V_m(expModes.good_ch{i},i);
    idx = find(V_eG == 1);
    V_mG = V_opt(expModes.good_ch{i},i);
    V_mG = V_mG / V_mG( idx);
    eigvecDiffOpt = [eigvecDiffOpt;(V_eG - V_mG) .* expModes.psiWeights(expModes.good_ch{i}, i)];
end

eigvecDiffOpt = eigvecDiffOpt(~isnan(eigvecDiffOpt));
f_opt = sqrt(Lambda_opt) / 2/ pi;




