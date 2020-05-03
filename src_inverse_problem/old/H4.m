function res = H4(l1, l2, ii, jj)
% Forth type of H function
% As parameters takes:
%   l1 - one of the boundaries (1 or 2)
%   l2 - one of the boundaries (1 or 2)
%   l1 != l2
%   t1, t2 - are from [0, 2pi]

%if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
%    error(['WRONG INDEX IN H4(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
%elseif l1 == l2
%    error(['l1 CANNOT BE EQUAL TO l2 IN H4(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
%end

global s;
global n_11;
global n_21;
global n_12;
global n_22;
global x1;
global y1;
global x2;
global y2; 

if l1 == 1 && l2 == 1
    n_1c = n_11(ii);
    n_2c = n_21(ii);
    
    n_1i = n_11(jj);
    n_2i = n_21(jj);
    
    x_c = x1(ii);
    y_c = y1(ii);
    
    x_i = x1(jj);
    y_i = y1(jj);
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
end

res = - (n_1c * (x_c - x_i) + n_2c * (y_c - y_i)) * (n_1i * (x_c - x_i) + n_2i * (y_c - y_i)) /...
         2 / r2(l1, l2, ii, jj)^2 - (n_1c * n_1i + n_2c * n_2i) * (1 + 2 * log(r2(l1, l2, ii, jj))) / 4;
end
