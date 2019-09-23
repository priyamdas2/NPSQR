function sum = log_like(Basis,linear_basis, theta,x_array,y_array)
m = size(theta,1) - 1;
p = size(theta,2) - 2;
n = max(size(x_array,1), size(x_array,2));
tau_seq = solve_tau(Basis,theta,x_array,y_array);
sum = 0;
tau_cum_probs = linspace(0,1,m+1);
differentiation = zeros(n,1);
for iii = 1:n
    %%% x - basis coeff calculation %%%%%%%%%%%%%
    x_value = x_array(iii);
    x_cum_probs = linspace(0,1,p+1);
    for j = 2:(p+1)
        if(x_value <= x_cum_probs(j))
            x_interval = j-1;
            break
        end
    end
    
    Basis_now_x = Basis(:,(3*(x_interval-1)+1) : (3*(x_interval-1)+3));
    Basis_values = Basis_now_x(:,1)*x_value^2 + Basis_now_x(:,2)*x_value + Basis_now_x(:,3);
    i_coeffs = theta*Basis_values; % delta i coeffs actually
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    for i = 2:(m+1)
        if(tau_seq(iii) <= tau_cum_probs(i)+0.000000001)
            tau_interval = i-1;
            break
        end
    end
    linear_basis_here = linear_basis(:,(2*(tau_interval-1)+1) : (2*(tau_interval-1)+2));
    differentiation(iii) = (m+2)*(transpose(i_coeffs)*linear_basis_here(:,1)*tau_seq(iii) + transpose(i_coeffs)*linear_basis_here(:,2));
    sum = sum - log(differentiation(iii));
end
end
