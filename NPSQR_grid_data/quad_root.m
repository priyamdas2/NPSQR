function [r_1,r_2]= quad_root(c,b,a)
if(a == 0)
    r_1 = -c/b;
    r_2 = -999;
else
    D = sqrt(b^2-4*a*c);
    r_1 = (-b-D)/(2*a);
    r_2 = (-b+D)/(2*a);
end
end