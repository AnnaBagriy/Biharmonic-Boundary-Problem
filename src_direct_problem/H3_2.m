function res = H3_2(l, ii, jj)
% Third type of H(2) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector
%   jj - index for s-vector
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H3_2(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

global s;
global n_11;
global n_21;
global n_12;
global n_22;
global x1;
global y1;
global x2;
global y2;

if ii == jj
    res = 0;
else
    if l == 1
            res = (n_11(ii) * (x1(ii) - x1(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((s(ii) - s(jj)) / 2))^2)) + 1)... 
                 + n_21(ii) * (y1(ii) - y1(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((s(ii) - s(jj)) / 2))^2)) + 1)) / 4;
    else
            res = (n_12(ii) * (x2(ii) - x2(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((s(ii) - s(jj)) / 2))^2)) + 1)... 
                 + n_22(ii) * (y2(ii) - y2(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((s(ii) - s(jj)) / 2))^2)) + 1)) / 4;
    end
end

end