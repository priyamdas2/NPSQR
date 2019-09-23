function answer = solve_tau(Basis,theta,x_array,y_array) % x_array = 0.7893; y_array =  0.5688
m = size(theta,1) - 1;
p = size(theta,2) - 2;
n = max(size(x_array,1), size(x_array,2));
answer = zeros(n,1);

for iii = 1:n
    x_value = x_array(iii);
    y_value = y_array(iii);
    
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
    
    tau_cum_probs = linspace(0,1,p+1);
    tau_critical_cumsum = zeros(m+1,1);
    tau_critical_cumsum(1) = 0;
    for i = 2:(m+1)
        Basis_now_tau = Basis(:,(3*(i-2)+1) : (3*(i-2)+3));
        tau_critical_cumsum(i) = [0,transpose(cumsum(i_coeffs))]*(Basis_now_tau(:,1)*tau_cum_probs(i)^2 + Basis_now_tau(:,2)*tau_cum_probs(i) + Basis_now_tau(:,3));
    end
    
    for i = 2:(m+1)
        if(y_value <= tau_critical_cumsum(i)+0.00000001)
            y_interval = i-1;
            break
        end
    end
    
    Basis_now_tau_again = Basis(:,(3*(y_interval-1)+1) : (3*(y_interval-1)+3));
    a = [0,transpose(cumsum(i_coeffs))]*Basis_now_tau_again(:,1); % coeff of tau^2
    a = round(a*10^4)/10^4;
    b = [0,transpose(cumsum(i_coeffs))]*Basis_now_tau_again(:,2); % coeff of tau
    c = [0,transpose(cumsum(i_coeffs))]*Basis_now_tau_again(:,3); % coeff of 1
    
    [r_1,r_2] = quad_root(c-y_value,b,a);
    if(r_1<=x_cum_probs(y_interval+1) && r_1>x_cum_probs(y_interval))
        answer(iii) = real(r_1);
    else
        answer(iii) = real(r_2);
    end
end


end




% m=5;
% p=5;
% starting_point = ones(m+1,p+2);
% for j = 1:(p+2)
%     starting_point(:,j) = starting_point(:,j) / sum(starting_point(:,j));
% end
% theta = starting_point;