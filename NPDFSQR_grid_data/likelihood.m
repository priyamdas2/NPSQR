function total_log_like = likelihood(c,Basis,coef_diff,x_array_transformed,positions,division_transformed)
data_length = max(size(x_array_transformed));
max_position = max(positions);

prob_division = zeros(data_length,1);

for i = 1:data_length
    if(positions(i) == 1)
        prob_division(i) = find_tau_eco(c,Basis,coef_diff,x_array_transformed(i), division_transformed(2));
    else
        if(positions(i) == max_position)
            prob_division(i) = 1 - find_tau_eco(c,Basis,coef_diff,x_array_transformed(i), division_transformed(max_position));
        else
            prob_division(i) = find_tau_eco(c,Basis,coef_diff,x_array_transformed(i), division_transformed(positions(i)+1))-...
                find_tau_eco(c,Basis,coef_diff,x_array_transformed(i), division_transformed(positions(i)));
        end
    end
end

total_log_like = sum(log(prob_division));
end
