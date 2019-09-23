function q_tau_x = q_tau_x_axis(c,Basis,theta,tau, x_axis)
p = size(c,1)-2;
m = size(theta,1)-1;

coef = vertcat(zeros(1,p+2),cumsum(theta));
for j = 1:m
    if(tau<=j/m)
        tau_loc = j;
        break
    end
end

for k = 1:p
    if(x_axis<=k/p)
        x_loc = k;
        break
    end
end

Basis_here = Basis(:,3*(tau_loc-1)+1 : 3*(tau_loc-1)+3);
c_here = c(:,3*(x_loc-1)+1 : 3*(x_loc-1)+3);

q_tau_x = 0;

for j=1:(m+2)
    for k = 1:(p+2)
        q_tau_x = q_tau_x + coef(j,k)*(Basis_here(j,1)*tau^2+Basis_here(j,2)*tau+Basis_here(j,3))*(c_here(k,1)*x_axis^2+c_here(k,2)*x_axis+c_here(k,3));
    end
end
    
