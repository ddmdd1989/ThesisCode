
% lsqnonlin
         
color = {'k','c','y','w'} ;
k = [-0.10 0.10] ;
figHand = figure; set (figHand, 'Position',[250 250 500 250]);
hold on
for j = 1 : 6
    for i = 1 :2
        h = bar(j + k(i),(error(i,j)),0.15,color{i}) ;
    end
end

grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itk}_1','{\itk}_2','{\itk}_3','{\itk}_4','{\itk}_5','{\itk}_6'} ;
h = my_xticklabels(gca,1:6,xtl);

lgd = legend('Case 2(a)','Case 2(b)');
set(gca,'FontSize',FZ);
set(lgd,'FontSize',FZ);

text(0.2,9e-13,'Case 2(a), Avg. Err. = 1.54\times10^{-12}%','Fontsize',FZ)
text(0.2,8e-13,'Case 2(b), Avg. Err. = 2.93\times10^{-12}%','Fontsize',FZ)
ylim([0,1e-12])




