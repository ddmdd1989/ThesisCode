function [M_C, K_C, P_K] = StrucCondens(M, K, P_K, K_al, M_al, M_C, n_s, n_i, n_r, n_q, condensationType)

% condensationType   1: Craig-Bampton condensation   else : Real-mode Consensation

N_tilde = size(P_K,1) ;
n_k = size(P_K,3) ;
N = size(M,1) ;

if condensationType == 1

    T = -K(n_s+n_i+1:N,n_s+n_i+1:N)\K(n_s+n_i+1:N,n_s+1:n_s+n_i) ;
    [Vs_raw,Ss_raw] = eig(K(n_s+n_i+1:N,n_s+n_i+1:N), M(n_s+n_i+1:N,n_s+n_i+1:N) ) ;
    [Ss,I] = sort((diag(Ss_raw)),'ascend') ;
    Vs = Vs_raw(:,I);                          % fix boundary mode shapes for residue structure

    H = [   eye(n_i)    zeros(n_i,n_q)
            T           Vs(:,1:n_q) ] ;
else
    % original strucuture modes
    [V_a,Lambda_a] = eig(K, M ) ;
    [Lambda_a,I] = sort((diag(Lambda_a)),'ascend') ;
    V_a = V_a(:,I);
    for i = 1 : n_q
        V_a(:,i) = V_a(:,i)/max(abs(V_a(:,i))) ;
    end
    
    H = zeros(n_i + n_r, n_i + n_q);
    H(1 : n_i, 1 : n_i) = eye(n_i);
    H(n_i+1 : n_i + n_r, n_i+1 :n_i + n_q) = V_a(n_s+n_i+1 : N, 1:n_q);
end


KR  = zeros(N_tilde) ;   
MR  = zeros(N_tilde) ;       
% CR  = zeros(n_s+ n_q) ;   

KR(n_s+1:N_tilde,n_s+1:N_tilde) = H'*K_al*H ;
MR(n_s+1:N_tilde,n_s+1:N_tilde) = H'*M_al*H ;
% CR(n_s:n_s+n_q,n_s:n_s+n_q) = M_T'*C_al*M_T ;


[VR_raw,SR_raw] = eig(KR(n_s+1:N_tilde,n_s+1:N_tilde),MR(n_s+1:N_tilde,n_s+1:N_tilde) ) ;
[SR,I] = sort((diag(SR_raw)),'ascend') ;
VR = VR_raw(:,I);
KR_diag = VR'*KR(n_s+1:N_tilde,n_s+1:N_tilde)*VR ;
VR_l = inv(VR) ;

n_zeros = 0 ;
for i = 1 : n_i+n_q
    if (SR(i)<1e-3)    
        n_zeros = n_zeros + 1 ;
    else
        break ;
    end
end

for i = n_zeros+1 : n_i + n_q
    P_K(n_s+1:N_tilde,n_s+1:N_tilde,n_k+i-n_zeros) = KR_diag(i,i)*VR_l(i,:)'*VR_l(i,:) ;
end

n_alpha = size(P_K,3);

% Condensed mass and stiffness matrics M_C and K_C

M_C = M_C+ MR ;
K_C = sum(P_K,3) ;
