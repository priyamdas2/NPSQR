function returned = fun(c,Basis,theta,x_array_transformed,positions,division_transformed) % theta has (m+1) rows and (p+2) columns
p= size(theta,2)-2;
coef_diff = cell(1,p+2);
for k = 1:(p+2)
    coef_diff{k} = theta(:,k);
end
returned = -likelihood(c,Basis,coef_diff,x_array_transformed,positions,division_transformed);
end
