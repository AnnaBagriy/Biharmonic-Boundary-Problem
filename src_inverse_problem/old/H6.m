function res = H6(l1, l2, ii, jj, nu)
% Sixth type of H function
% As parameters takes:
%   l1 - one of the boundaries (1 or 2)
%   l2 - one of the boundaries (1 or 2)
%   t1, t2 - are from [0, 2pi]

global n_11;
global n_21;
global n_12;
global n_22;

global x1;
global y1;
global x2;
global y2;


global x_derivative1;
global x_derivative2;
global y_derivative1;
global y_derivative2;

global x_2derivative1;
global x_2derivative2;
global y_2derivative1;
global y_2derivative2;

%TODO: Refactor

if l1 == 1 && l2 == 1
    n_1c = n_11(ii);
    n_2c = n_21(ii);
    
    n_1i = n_11(jj);
    n_2i = n_21(jj);
    
    x_c = x1(ii);
    y_c = y1(ii);
    
    x_i = x1(jj);
    y_i = y1(jj);
    
    x_1der_1 = x_derivative1(ii);
    x_1der_2 = y_derivative1(ii);

    x_2der_1 = x_2derivative1(ii);
    x_2der_2 = y_2derivative1(ii);
end
if l1 == 1 && l2 == 2
    n_1c = n_11(ii);
    n_2c = n_21(ii);
    
    n_1i = n_12(jj);
    n_2i = n_22(jj);
    
    x_c = x1(ii);
    y_c = y1(ii);
    
    x_i = x2(jj);
    y_i = y2(jj);
    
    x_1der_1 = x_derivative1(ii);
    x_1der_2 = y_derivative1(ii);

    x_2der_1 = x_2derivative1(ii);
    x_2der_2 = y_2derivative1(ii);
end
if l1 == 2 && l2 == 1
    n_1c = n_12(ii);
    n_2c = n_22(ii);
    
    n_1i = n_11(jj);
    n_2i = n_21(jj);
    
    x_c = x2(ii);
    y_c = y2(ii);
    
    x_i = x1(jj);
    y_i = y1(jj);
    
    x_1der_1 = x_derivative2(ii);
    x_1der_2 = y_derivative2(ii);

    x_2der_1 = x_2derivative2(ii);
    x_2der_2 = y_2derivative2(ii);
end
if l1 == 2 && l2 == 2
    n_1c = n_12(ii);
    n_2c = n_22(ii);
    
    n_1i = n_12(jj);
    n_2i = n_22(jj);
    
    x_c = x2(ii);
    y_c = y2(ii);
    
    x_i = x2(jj);
    y_i = y2(jj);
    
    x_1der_1 = x_derivative2(ii);
    x_1der_2 = y_derivative2(ii);

    x_2der_1 = x_2derivative2(ii);
    x_2der_2 = y_2derivative2(ii);
end

if ii == jj
    res = (1 - 3 * nu) * (n_1c * x_2der_1 + n_2c * x_2der_2) /...
        (4 * (x_1der_1^2 + x_1der_2^2));
else
    res = (1 - nu) * (((n_1c * (x_c - x_i) + n_2c *...
        (y_c - y_i))^2 * (n_1i * (x_c - x_i) +...
        n_2i * (y_c - y_i))) / r2(l1, l2, ii, jj)^4 -...
        ((n_1c * n_1i + n_2c * n_2i) * (n_1c *...
        (x_c - x_i) + n_2c * (y_c - y_i))) /...
        r2(l1, l2, ii, jj)^2) / (2 * pi) - ((1 + nu) * (n_1i *...
        (x_c - x_i) + n_2i * (y_c - y_i)) /...
        (4 * pi * r2(l1, l2, ii, jj)^2));
end

end