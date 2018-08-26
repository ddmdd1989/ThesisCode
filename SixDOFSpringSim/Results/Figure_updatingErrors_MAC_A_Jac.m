
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
FZ = 12;
grid on

ylabel('Relative error \fontname{Times New Roman}{\ite_i}(%)','fontsize',FZ)
xtl = {'{\itk}_1','{\itk}_2','{\itk}_3','{\itk}_4','{\itk}_5','{\itk}_6'} ;
h = my_xticklabels(gca,1:6,xtl);

lgd = legend('Case 1(a)','Case 1(b)');
set(gca,'FontSize',FZ);
set(lgd,'FontSize',FZ);
text(0.2,2.8e-5,'Case 1(a), Avg. Err. = 2.34\times10^{-6}%','Fontsize',FZ)
text(0.2,2.5e-5,'Case 1(b), Avg. Err. = 9.74\times10^{-6}%','Fontsize',FZ)
ylim([0,3e-5])



