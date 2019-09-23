function answer = likelihood_dist(Basis, theta,x_array_transformed,positions,division_transformed)
data_length = max(size(x_array_transformed));
max_position = max(positions);
log_sum_tot = zeros(data_length,1);
for iii = 1:data_length
    x_value = x_array_transformed(iii);
    if(positions(iii) == 1)
        y_value = division_transformed(2);
        log_sum_tot(iii) = log(like_dist_DIST(Basis, theta, x_value, y_value));
    else
        if(positions(iii) == max_position)
            y_value = division_transformed(max_position);
            log_sum_tot(iii) =  log(1-like_dist_DIST(Basis, theta, x_value, y_value));
        else
            y_value_2 = division_transformed(positions(iii)+1);
            y_value_1 = division_transformed(positions(iii));
            log_sum_tot(iii) = log(like_dist_DIST(Basis, theta, x_value, y_value_2) - ...
                like_dist_DIST(Basis, theta, x_value, y_value_1));
        end
    end
end
answer = sum(log_sum_tot);
end
% x_array = 0.45;
% y_array = 0.89;