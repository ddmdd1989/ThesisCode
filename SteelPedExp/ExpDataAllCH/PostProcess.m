%% 
clear
close all
LoadStructure;

n_alpha = size(K_j,3) - 2;
Modeindex = [1 2 4];

%% Assemble structure parameter:
structModel.M0 = M0;
structModel.K0 = K0;
structModel.K_j = K_j(:,:,1:n_alpha);


%% Experimental modal properties
load ModeInfoAllCh

badCh = setdiff(1:length(measDOFs),good_ch{1});
measDOFs(badCh) = [];
mode_mean(badCh,:) = [];

numModes = length(Modeindex);


% 
Lambda_e = Lambda_mean(Modeindex);
UDOF = setdiff(1:N,measDOFs);
reord = [measDOFs';UDOF'];
K0 = K0(reord,reord);
M0 = M0(reord,reord);
K_j = K_j(reord,reord,:);
[V,D] = eigs(K0,M0,20,'sm');
[Lambda_a,srt_idx] = sort(diag(D),'ascend');
V_ma = V(1:length(measDOFs),srt_idx);
V_ma = V_ma(:,Modeindex);
Lambda_a = Lambda_a(Modeindex);
f_e = sqrt(Lambda_e) / 2 / pi;
f_a = sqrt(Lambda_a) / 2/ pi;

MAC_ini = mac(V_ma,mode_mean(:,Modeindex));
load Space_form2_JACon_TRR


[fval_sort,index_sort] = sort(fval);

figHand = figure; 

set (figHand, 'Position',[250 250 500 250]);
semilogy(index_sort(1:52),fval(index_sort(1:52)),'*k')
hold on 
semilogy(index_sort(53:end),fval(index_sort(53:end)) ,'ok')
xlabel('Starting Point #','Fontsize',11);
ylabel('Obj. Func. Value','Fontsize',11)
hold off

% Plot alpha value
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : n_alpha    
    h = bar(j,alpha(j,index_sort(1)),0.15,'k');
end

FZ = 10.5;
grid on
ylabel('Relative change \fontname{Times New Roman}{\it\alpha_i}','fontsize',FZ)
xtl = {'{\it\alpha}_1','{\it\alpha}_2','{\it\alpha}_3','{\it\alpha}_4','{\it\alpha}_5','{\it\alpha}_6',...
       '{\it\alpha}_{7}','{\it\alpha}_{8}','{\it\alpha}_{9}','{\it\alpha}_{10}','{\it\alpha}_{11}','{\it\alpha}_{12}'...
       '{\it\alpha}_{13}','{\it\alpha}_{14}'};
h = my_xticklabels(gca,1:n_alpha,xtl);
text(0.2,0.8,'Obj. Val = 1.3117\times10^3','fontsize',FZ)[~,idx] = min(fval);
opt = alpha(:,idx);


K = sum(K_j,3);
for i = 1:n_alpha
    K = K + K_j(:,:,i) * opt(i);
end

[V,D] = eigs(K,M0,20,'sm');
[Lambda_opt,srt_idx] = sort(diag(D),'ascend');
V_opt = V(1:length(measDOFs),srt_idx);
V_opt = V_opt(:,Modeindex);
Lambda_opt = Lambda_opt(Modeindex);
f_opt = sqrt(Lambda_opt) / 2/ pi;
MAC_opt = mac(mode_mean(:,Modeindex),V_opt);



ylim([-1,1]);
hold off

