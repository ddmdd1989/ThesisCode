clc;clear;close all

alpha_act = [-0.25; -0.20; -0.15; -0.10; -0.05; 0.05];
         
n_alpha = length(alpha_act);

num_star = 100;

%% MAC Jacon
filename = 'LumpedSpringMass_form1_JACon_LM';
load(filename);
[fval_sort1,index] = sort(fval,'ascend');
x = alpha;
x_sort = alpha(:,index);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;

for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end


filename = 'LumpedSpringMass_form1_JACon_TRR';
load(filename);
[fval_sort2,index] = sort(fval,'ascend');
x = alpha;
x_sort = alpha(:,index);
error(2,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;

for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,2) = mean(error_iter);
end



Figure_updatingErrors_MAC_A_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite_e}_{avg}(%)','fontsize',FZ);
lgd = legend('Case 1(a)','Case 1(b)');
set(gca,'FontSize',FZ);
set(lgd,'FontSize',FZ);




%% Eigval Diff Ana Jac
filename = 'LumpedSpringMass_form2_JACon_LM';
load(filename);
[fval_sort5,index] = sort(fval,'ascend');
x = alpha;
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
 
for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end 


filename = 'LumpedSpringMass_form2_JACon_TRR';
load(filename);
x = alpha;
[fval_sort6,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(2,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;


for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,2) = mean(error_iter);
end

Figure_updatingErrors_Vm_A_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite_e}_{avg}(%)','fontsize',FZ);
lgd = legend('Case 2(a)','Case 2(b)');
set(gca,'FontSize',FZ);
set(lgd,'FontSize',FZ);





%% Modal Dynamic Residual Ana Jac
filename = 'LumpedSpringMass_form3_JACon_LM';
load(filename);
x = alpha(1:n_alpha,:);
[fval_sort7,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(1,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
 
for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,1) = mean(error_iter);
end 


filename = 'LumpedSpringMass_form3_JACon_TRR';
load(filename);
x = alpha(1:n_alpha,:);
[fval_sort8,index] = sort(fval,'ascend');
x_sort = x(:,index);
n_alpha = size(x,1);
error(2,:) = abs(x_sort(:,1) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;


for i = 1:num_star
    error_iter = abs(x(:,i) - alpha_act) ./ (ones(n_alpha,1) + alpha_act) * 100;
    mean_error(i,2) = mean(error_iter);
end

Figure_updatingErrors_Res_A_Jac


figHand = figure; 
set (figHand, 'Position',[250 250 500 250]);

plot(1:num_star,mean_error(:,1),'*','Color', [17/255 17/255 17/255])
hold on
plot(1:num_star,mean_error(:,2),'ok')
xlabel('Starting Points','Fontsize',11);
ylabel('Avg. relative error \fontname{Times New Roman}{\ite_e}_{avg}(%)','fontsize',FZ);
lgd = legend('Case 3(a)','Case 3(b)');
set(gca,'FontSize',FZ);
set(lgd,'FontSize',FZ);













      