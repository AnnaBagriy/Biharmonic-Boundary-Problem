function res = H4_1(l, ii, jj)
% Fourth type of H(1) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector
%   jj - index for s-vector
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

global n_11;
global n_21;
global n_12;
global n_22;

if l == 1
    res = - (n_11(ii) * n_11(jj) + n_21(ii) * n_21(jj)) / 4;
else
    res = - (n_12(ii) * n_12(jj) + n_22(ii) * n_22(jj)) / 4;  
end


end