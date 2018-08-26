function flag = ConvergenceCheck(obj_old, obj_new, x_inc, x_new, TolFun, TolX)

% function flag = ConvergenceCheck(obj_old, obj_new, x_inc, x_new, TolFun, TolX)
%
%	(c) Yang Wang, Daepng Zhu, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function determines the termination of the optimization process based on the 
% change of objective function and optimization variables values.
%
% Input:
%   obj_old: objective function value of last iteration
%	obj_new: objective function value of current iteration
%	x_inc: change of optimization variables values from last iteration to
%		   current iteration
%	x_new: optimization varables vaules of current iteration
%   TolX: Termination tolerance on the value change of the optimization
%   	  variables between two iterations; 
%   TolFun: Termination tolerance on the value change of the objective
%      	    function between two iterations;
% Output:
% 	flag: 0 - continue next iteration of the optimization process
%		  1 - terminate the optimization process

diffObj = abs(obj_new-obj_old) ;
diffX = norm(x_inc) ;
  
if (diffObj<=TolFun * (1+abs(obj_old)) || diffX<=TolX *(1+norm(x_new)) )
      flag = 1 ;
else
      flag = 0;
end
    