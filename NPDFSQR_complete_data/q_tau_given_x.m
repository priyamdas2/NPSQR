function value = q_tau_given_x(Basis, theta,tau_value,x_value)
m = size(theta,1) - 1;
p = size(theta,2) - 2;

x_cum_probs = linspace(0,1,p+1);
for j = 2:(p+1)
    if(x_value <= x_cum_probs(j))
        x_interval = j-1;
        break
    end
end

Basis_now_x = Basis(:,(3*(x_interval-1)+1) : (3*(x_interval-1)+3));
Basis_values = Basis_now_x(:,1)*x_value^2 + Basis_now_x(:,2)*x_value + Basis_now_x(:,3);

i_coeffs = theta*Basis_values;

tau_cum_probs = linspace(0,1,m+1);
for i = 2:(p+1)
    if(tau_value <= tau_cum_probs(i)+0.0000001)
        tau_interval = i-1;
        break
    end
end

Basis_now_tau_again = Basis(:,(3*(tau_interval-1)+1) : (3*(tau_interval-1)+3));
cumsum_i_coeffs = [0,transpose(cumsum(i_coeffs))];

value = cumsum_i_coeffs*Basis_now_tau_again(:,1)*tau_value^2 + ...
    cumsum_i_coeffs*Basis_now_tau_again(:,2)*tau_value + ...
    cumsum_i_coeffs*Basis_now_tau_again(:,3);
end
