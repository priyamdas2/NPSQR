function e = inv_pareto_dist(z)
a = 0.45;
k = 4.9;
s = 52;
n = size(z,2);
e = zeros(1,n);
for i=1:n
    e(i) = s*((1/(1-z(i))^(1/a))-1)^(1/k);
end
end

