function res = H2_1(l, ii, jj)
% Second type of H(1) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector 
%   jj - index for s-vector 
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H2_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

global n_12;
global n_22;
global n_11;
global n_21;
global x1;
global y1;
global x2;
global y2;

if ii == jj
    res = 0;
else
    if l == 1
        res = -(n_11(jj) * (x1(ii) - x1(jj)) + n_21(jj) * (y1(ii) - y1(jj))) / 4;
    else
        res = -(n_12(jj) * (x2(ii) - x2(jj)) + n_22(jj) * (y2(ii) - y2(jj))) / 4;
    end
end

end