function res = H5_1(l, ii, jj, nu)
% Fifth type of H(1) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector
%   jj - index for s-vector
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

res = (1 + nu) / 4;

end