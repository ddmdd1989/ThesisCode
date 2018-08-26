
% lsqnonlin
         
color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 15
    for i = 1 :2
        h = bar(j + k(i),(error(i,j)),0.15,color{i}) ;
    end
end

grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itE}_1','{\itE}_2','{\itE}_3','{\itE}_4','{\itE}_5','{\itE}_6',...
       '{\itE}_{t2}','{\itE}_{t3}','{\itE}_{t4}','{\itE}_{t5}','{\itE}_{t6}',...
       '{\itk}_{y1}','{\itk}_{z1}','{\itk}_{y2}','{\itk}_{z2}'} ;
h = my_xticklabels(gca,1:15,xtl);

lgd = legend('Case 3(a)','Case 3(b)');
set(lgd,'fontsize',FZ);
set(gca,'fontsize',FZ);
text(0.2,2.7e-11,'Case 3(a), Avg. Err. = 3.11\times10^{-12}%','Fontsize',FZ)
text(0.2,2.4e-11,'Case 3(b), Avg. Err. = 5.17\times10^{-12}%','Fontsize',FZ)
ylim([0,3e-11])
% 


