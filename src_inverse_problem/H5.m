function res = H5(l1, l2, ii, jj)
% Fifth type of H function
% As parameters takes:
%   l1 - one of the boundaries (1 or 2)
%   l2 - one of the boundaries (1 or 2)
%   l1 != l2
%   t1, t2 - are from [0, 2pi]

global nu;

global n_11;
global n_21;
global n_12;
global n_22;

global x1;
global y1;
global x2;
global y2;

%TODO: Refactor

if l1 == 1 && l2 == 1
    n_1c = n_11(ii);
    n_2c = n_21(ii);
    
    x_c = x1(ii);
    y_c = y1(ii);
    
    x_i = x1(jj);
    y_i = y1(jj);
end
if l1 == 1 && l2 == 2
    n_1c = n_11(ii);
    n_2c = n_21(ii);
    
    x_c = x1(ii);
    y_c = y1(ii);
    
    x_i = x2(jj);
    y_i = y2(jj);
end
if l1 == 2 && l2 == 1
    n_1c = n_12(ii);
    n_2c = n_22(ii);
   
    x_c = x2(ii);
    y_c = y2(ii);
    
    x_i = x1(jj);
    y_i = y1(jj);
end
if l1 == 2 && l2 == 2
    n_1c = n_12(ii);
    n_2c = n_22(ii);
    
    x_c = x2(ii);
    y_c = y2(ii);
    
    x_i = x2(jj);
    y_i = y2(jj);
end

res = (1 + 3 * nu) / 8 * pi + ((1 - nu) * (n_1c * (x_c - x_i) +...
    n_2c * (y_c - y_i))^2) / ( 4 * pi * r2(l1, l2, ii, jj)^2) +...
    ((1 + nu) * log(r2(l1, l2, ii, jj)^2)) / (8 * pi);

end

