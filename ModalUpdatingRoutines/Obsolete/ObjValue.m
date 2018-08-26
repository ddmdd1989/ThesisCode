function y = ObjValue(M,K,Lambda_e,V_ep,m,weight,t,alphasum_old)

for i = 1 : m
    tt(i) = norm((-M * Lambda_e(i)  + K )* V_ep(:,i))*weight(i);
end

y = norm(tt,2)^2+t^2*norm(alphasum_old,2)^2 ;