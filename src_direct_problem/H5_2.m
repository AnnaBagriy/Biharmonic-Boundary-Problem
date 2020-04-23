function res = H5_2(l, ii, jj)
% Fifth type of H(2) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector
%   jj - index for s-vector
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

global nu;

global s;

global x_derivative1;
global x_derivative2;
global y_derivative1;
global y_derivative2;

if l == 1
     if ii == jj
         res = (1 + 3 * nu) / 4 + (1 + nu) * log(exp(1) * (x_derivative1(ii)^2 + y_derivative1(ii)^2)^2) / 4;
     else
         res = H5(l, l, ii, jj) - (1 + nu) * log(4 * (sin((s(ii) - s(jj)) / 2))^2  / exp(1)) / 4;
     end
else
     if ii == jj
         res = (1 + 3 * nu) / 4 + (1 + nu) * log(exp(1) * (x_derivative2(ii)^2 + y_derivative2(ii)^2)^2) / 4;
     else
         res = H5(l, l, ii, jj) - (1 + nu) * log(4 * (sin((s(ii) - s(jj)) / 2))^2  / exp(1)) / 4;
     end
end

end