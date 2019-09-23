function answer = DIST_marginal_tau_given_x(Basis, theta, taus, x_value,size_grids_for_cdf) % x_value in [0,1]; x_value = 0.65

grids_for_cdf = linspace(0,1,size_grids_for_cdf);

check_values = zeros(size_grids_for_cdf,1);
for ii = 1:size_grids_for_cdf
    y_value = grids_for_cdf(ii);
    check_values(ii) = y_given_x_dist(Basis, theta, x_value, y_value);
end
tau_grid_size = max(size(taus));
corresponding_y_value = zeros(tau_grid_size,1);

for jj = 1:tau_grid_size
    for ii = 1:size_grids_for_cdf
        if(taus(jj)<=check_values(ii))
            corresponding_y_value(jj) = ((check_values(ii)-taus(jj))*(grids_for_cdf(ii-1)) + ...
                (taus(jj) - check_values(ii-1))*(grids_for_cdf(ii)))/(check_values(ii)-check_values(ii-1));
            break
        end
    end
end
answer = corresponding_y_value;
end