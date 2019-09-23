function answer = log_like_dist(Basis,linear_basis, theta, x_array, y_array)
m = size(theta,1) - 1; % y basis coeff diffs 
p = size(theta,2) - 2; % x basis coeffs
n = max(size(x_array,1), size(x_array,2));
log_sum_tot = zeros(n,1);
for iii = 1:n
    x_value = x_array(iii);
    y_value = y_array(iii);
    % basis value at y
    for i = 1:m
        if(y_value <= i/m)
            y_interval = i;
            break
        end
    end
    
    linear_basis_y_now = linear_basis(:,(2*(y_interval-1)+1): (2*(y_interval-1)+2));
    
    % basis value at x
    for j = 1:p
        if(x_value <= j/p)
            x_interval = j;
            break
        end
    end
    
    Basis_x_now = Basis(:,(3*(x_interval-1)+1): (3*(x_interval-1)+3));
    
    value = zeros(p+2,1);
    sum_tot = 0;
    for j = 1:(p+2)
        value(j) = (m+2)*(transpose(theta(:,j))*linear_basis_y_now(:,1)*y_value + ...
            transpose(theta(:,j))*linear_basis_y_now(:,2));
        B_j_x = sum((value(j)*Basis_x_now(j,:)).*[x_value^2, x_value, 1]);
        sum_tot = sum_tot+B_j_x;
    end
    log_sum_tot(iii) = log(sum_tot);
end

answer = sum(log_sum_tot);
end
% x_array = 0.45;
% y_array = 0.89;