function res = H1_2(l, ii, jj)
% First type of H(2) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector
%   jj - index for s-vector
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H1_2(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

global s;

if ii == jj
    res = 0;
else
    res = r2(l, l, ii, jj).^2 .* log((exp(1) .* r2(l, l, ii, jj)^2) ./ (4 .* (sin((s(ii) - s(jj)) ./ 2)).^2)) ./ 8;
end

end