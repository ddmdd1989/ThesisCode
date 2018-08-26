function HyperlineWalk(fun, xStart, xInput, numStep, ub, lb, opt)
%% function HyperlineWalk(fun, xStart, xIput2, num_step, ub, lb, opt)
%
%
%   (c) Yang Wang, Xinjun Dong, Dan Li (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
%   This function is to draw a hyperline walk from Xstart to the boundary
%   of updating variables. 
%%% Input:
%     fun - objective function
%     xStart - starting point of the hyperline walk
%     xInput:
%       if (opt = 1) - second fixed point on the hyperline walk
%       else         - gradient of the hyperline walk at starting point
%     num_step:
%       if (opt = 1) - number of steps between starting and second fixed
%                      point on the hyperline walk
%       else         - not used
%     ub: upper bound of updating variables
%     lb: lower bound of updating variables
%     opt: 1 - hyperline walk crossing two fixed points (xStart xInput)
%          2 - hyperline walk along a fixed deirection (xInput) from the
%              strating point(xStart)
xPoint = [];
t = 0;
if(opt == 1)
    grad = (xInput - xStart) / numStep;
else
    grad = xInput;
end
while(1)
    xPoint_back = xStart + t * grad;
    if(isempty( find(xPoint_back < lb, 1))&& isempty(find(xPoint_back > ub, 1)))
        xPoint = [xPoint_back xPoint];
        t = t - 1;
    else
        break;
    end
    
end
t_min = t;
t = 1;
while(1)
    xPoint_fw = xStart + t * grad;
    if(isempty( find(xPoint_fw < lb, 1))&& isempty(find(xPoint_fw > ub, 1)))
        xPoint = [xPoint xPoint_fw];
        t = t + 1;
    else
        break;
    end
    
end

t_sp = t_min+1:1:t-1;
for i = 1:size(xPoint,2)
   
   PDiff(i) = norm(feval(fun,xPoint(:,i)),2);
    
end

figure
semilogy(t_sp, PDiff,'k','LineWidth',1.5);
set(gca,'FontSize',10);
xlabel('Step','Fontsize',11);
ylabel('Obj. Func. Value','Fontsize',11);

end