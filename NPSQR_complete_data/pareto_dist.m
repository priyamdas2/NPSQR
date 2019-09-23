function e = pareto_dist(y)
a = 0.45;
k = 4.9;
s = 52;
n = size(y,1);
e = zeros(n,1);
for i = 1:n
    e(i) = 1-1/(1+(y(i)/s)^k)^a;
end
end





