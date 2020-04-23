function res = H2(l1, l2, ii, jj)
% Second type of H function
% As parameters takes:
%   l1 - one of the boundaries (1 or 2)
%   l2 - one of the boundaries (1 or 2)
%   l1 != l2
%   t1, t2 - are from [0, 2pi]

if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
    error(['WRONG INDEX IN H2(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
elseif l1 == l2
    error(['l1 CANNOT BE EQUAL TO l2 IN H2(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
end

global n_12;
global n_22;
global n_11;
global n_21;
global x1;
global y1;
global x2;
global y2;

if l2 == 1
    res = - (n_11(jj) * (x2(ii) - x1(jj)) * (1 + 2 * log(r2(l1, l2, ii, jj))) +...
             n_21(jj) * (y2(ii) - y1(jj)) * (1 + 2 * log(r2(l1, l2, ii, jj)))) / 4;
else
    res = - (n_12(jj) * (x1(ii) - x2(jj)) * (1 + 2 * log(r2(l1, l2, ii, jj))) +...
             n_22(jj) * (y1(ii) - y2(jj)) * (1 + 2 * log(r2(l1, l2, ii, jj)))) / 4;
end

end