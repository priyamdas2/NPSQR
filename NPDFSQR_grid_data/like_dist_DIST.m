function answer = like_dist_DIST(Basis, theta, x_value, y_value)
m = size(theta,1) - 1; % y basis coeff diffs 
p = size(theta,2) - 2; % x basis coeffs

    % basis value at y
    for i = 1:m
        if(y_value <= i/m)
            y_interval = i;
            break
        end
    end
    
    Basis_y_now = Basis(:,(3*(y_interval-1)+1): (3*(y_interval-1)+3));
    
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
        cumsum_theta_coeff = cumsum([0; theta(:,j)]);
        value(j) = (transpose(cumsum_theta_coeff)*Basis_y_now(:,1)*y_value^2 + ...
            transpose(cumsum_theta_coeff)*Basis_y_now(:,2)*y_value + ...
        transpose(cumsum_theta_coeff)*Basis_y_now(:,3));
        B_j_x = sum((value(j)*Basis_x_now(j,:)).*[x_value^2, x_value, 1]);
        sum_tot = sum_tot+B_j_x;
    end

answer = sum_tot;
end
% x_array = 0.45;
% y_array = 0.89;