function tau_value_array = find_tau_eco(c,Basis,coef_diff,x_axis, value_y_array) % x_axis =  0.5388; value_y_array = 0.3786
p = size(c,1)-2;
m = size(coef_diff{1},1)-1;
coeff = cell(p+2,1);
for k = 1:(p+2)
    coeff{k} = [0;cumsum(coef_diff{k})];
end

array = repmat(1/p,p,1);
cum_array = cumsum(array);

if(x_axis == 1)
    x_loc = p;
else
    for k=1:p
        if(x_axis<=cum_array(k))
            x_loc = k;
            break
        end
    end
end

c_basis_now = c(:,3*(x_loc-1)+1 : 3*(x_loc-1)+3);

border_given_x_axis = zeros(m,1);

for i = 1:m
    for j = 1:(m+2)
        for k = 1:(p+2)
            border_given_x_axis(i) = border_given_x_axis(i) + coeff{k}(j)*(c_basis_now(k,1)*x_axis^2 + c_basis_now(k,2)*x_axis +  c_basis_now(k,3))*...
                (Basis(j,3*(i-1)+1)*(i/m)^2 + Basis(j,3*(i-1)+2)*(i/m)+  Basis(j,3*(i-1)+3));
        end
    end
end

array_size = max(size(value_y_array,1),size(value_y_array,2));

tau_value_array = zeros(1,array_size);



for iii = 1:array_size
    value_y = value_y_array(iii);
    
    for ii = 1:m
        if (value_y <= border_given_x_axis(ii))
            value_y_location = ii;
            break
        end
    end
    
    Basis_now = Basis(:,3*(value_y_location-1)+1 : 3*(value_y_location-1)+3);
    
    coeff_of_x_square = 0;
    coeff_of_x = 0;
    coeff_of_1 = 0;
    for j = 1:(m+2)
        for k = 1:(p+2)
            coeff_of_x_square = coeff_of_x_square + coeff{k}(j)*(c_basis_now(k,1)*x_axis^2 + c_basis_now(k,2)*x_axis +  c_basis_now(k,3))*Basis_now(j,1);
            coeff_of_x = coeff_of_x + coeff{k}(j)*(c_basis_now(k,1)*x_axis^2 + c_basis_now(k,2)*x_axis +  c_basis_now(k,3))*Basis_now(j,2);
            coeff_of_1 = coeff_of_1 + coeff{k}(j)*(c_basis_now(k,1)*x_axis^2 + c_basis_now(k,2)*x_axis +  c_basis_now(k,3))*Basis_now(j,3);
        end
    end
    % Assuming equation as  a.x^2+b.x+c=0
    a = round(coeff_of_x_square*10000)/10000;
    b = round(coeff_of_x*10000)/10000;
    c = round((coeff_of_1-value_y)*10000)/10000;
    
    [r_1,r_2] = quad_root(c,b,a);
    
    %     if(r_1>0.99999 && r_1<1.0000001)
    %         ta = 1;
    %     else
    %         if(r_1 <= value_y_location/m && r_1 > (value_y_location-1)/m)
    %             ta = r_1;
    %         else
    %             ta = min(1,r_2);
    %         end
    %     end
    
    
    if(r_1 <= value_y_location/m+0.02 && r_1 > (value_y_location-1)/m-0.02)
        ta = r_1;
    else
        ta = min(1,r_2);
    end
    
    
    tau_value = real(ta);
    tau_value_array(iii) = tau_value;
end
end





