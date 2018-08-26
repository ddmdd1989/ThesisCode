function alphasum = ItertiveLsqnonlinModalMass(M_C, K_C, Lambda_e, V_m, P_K,P_M, weight, t, max_iter, TolFun, TolX, modalExpansionType,MDOF,UDOF)
 
% function alphasum = ItertiveLSQ_modalMass(M_C, K_C, Lambda_e, V_m, P, weight, t, max_iter, 
%                                       TolFun, TolX, modalExpansionType,MDOF,UDOF)
%
%   (c) Yang Wang, Dapeng Zhu, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2016
%
% This function tries to solve the optimization problem of modal dynamic residual approach 
% through two-step iterative optimization process
% Note: this function intends to update stiffness and modal mass of the
% structure; modal mass is only used in the substructure model updating
%
% Input:
%   M_C: mass matrix;
%   K_C: stiffness matrix;
%   Lambda_e: experimental eigenvalues;
%   V_m: experimental eigenvectors at measured DOFs;
%   P_K: influence matrices of stiffness variables;
%   P_M: influence matrices of modal mass variables;
%   weight: weighing factor for each experimentally-measured modes;
%   t: Laplace multiplier to avoid singular matrix inversion;
%   max_iter: maximum number of iterations;
%   TolX: Termination tolerance on the value change of the optimization
%         variables between two iterations; 
%   TolFun: Termination tolerance on the value change of the objective
%           function between two iterations;
%   modalExpansionType:
%           1: Full dynamic expansion       
%           2: Partial dynamic expansion
%   MDOF: Measured DOFs
%   UDOF: Unmeasured DOFs 
% Output:
%   alphasum: optimal values of structural parameters




N_tilde = size(M_C,1) ;
n_alpha = size(P_K,3) ;
n_beta = size (P_M,3) ;
numModes = size(V_m,2) ;
V_ep = zeros(N_tilde,numModes) ;
alphasum = zeros(n_alpha+n_beta,1) ;  


for ii = 1 : max_iter
    
    if (modalExpansionType == 1)  % least square expansion    
        for i = 1 : numModes
            Ai = [     (-Lambda_e(i)*M_C(MDOF,UDOF)+K_C(MDOF,UDOF))
                (-Lambda_e(i)*M_C(UDOF,UDOF)+K_C(UDOF,UDOF))    ];
            bi = [     (-Lambda_e(i)*M_C(MDOF,MDOF)+K_C(MDOF,MDOF))
                (-Lambda_e(i)*M_C(UDOF,MDOF)+K_C(UDOF,MDOF))    ];
            CT = -Ai\bi ;
            V_ep(MDOF,i) = V_m(:,i) ;
            V_ep(UDOF,i) = CT * V_m(:,i) ;
        end
    else
        % dynamic expansion
         for i = 1 : numModes
             CT = -(-Lambda_e(i)*M_C(UDOF,UDOF)+K_C(UDOF,UDOF))\(-Lambda_e(i)*M_C(UDOF,MDOF)+K_C(UDOF,MDOF)) ;
             V_ep(MDOF,i) = V_m(:,i) ;
             V_ep(UDOF,i) = CT * V_m(:,i) ;
         end
    end
    
    if ii ==1
        alphasum_old = zeros(n_alpha,1) ;
        obj_old = ObjValue(M_C,K_C,Lambda_e,V_ep,numModes,weight,t,alphasum_old) ;       
    else
        obj_old = obj_new ;
        alphasum_old = alphasum ;
    end
    
    % LSQ
    for j = 1 : numModes
        for i = 1 : n_alpha
            B((j-1)*N_tilde+1:j*N_tilde, i) = P_K(:,:,i) * V_ep(:,j) * weight(j) ;
        end
        for i = 1 : n_beta
            B((j-1)*N_tilde+1:j*N_tilde, i+n_alpha) = -Lambda_e(j) * P_M(:,:,i) * V_ep(:,j) * weight(j) ;
        end
        b((j-1)*N_tilde+1:j*N_tilde,1) = -(K_C-Lambda_e(j)*M_C)*V_ep(:,j) * weight(j);
    end
    
    alphaInc = (B'*B+t^2*eye(n_alpha+n_beta))\ (B'*b) ;
    alphasum = alphasum + alphaInc ;
    
    for i = 1 : n_alpha
        K_C = K_C + alphaInc(i) * P_K(:,:,i) ;
    end
    for i = 1 : n_beta
        M_C = M_C + alphaInc(i+n_alpha) * P_M(:,:,i) ;
    end
    
    
    obj_new = ObjValue(M_C,K_C,Lambda_e,V_ep,numModes,weight,t,alphasum) ;
 
    if (ConvergenceCheck(obj_old, obj_new, alphaInc, alphasum, TolFun, TolX))
          break
    end
    
    if ii == max_iter
            disp('Warning: Exceed the maximum iteration')
    end
end   
