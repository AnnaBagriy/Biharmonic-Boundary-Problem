function res = H4_2(l, ii, jj)
% Fourth type of H(2) function
% As parameters takes:
%   l - boundaries (1 or 2)
%   l1 = l2 = l
%   ii - index for s-vector
%   jj - index for s-vector
%   s - vector

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
end

global s;
%global n_11;
%global n_21;
%global n_12;
%global n_22;
%global x1;
%global y1;
%global x2;
%global y2;

if ii == jj
    syms t;
    y_derivative = matlabFunction(diff(y_vector(l,t),t));
    x_derivative = matlabFunction(diff(x_vector(l,t),t));
    res = -(1+log(exp(1)*((y_derivative(s(ii)))^2+(x_derivative(s(ii)))^2)))/4;
else
    res = H4(l,l, ii, jj) - H4_1(l, ii, jj)* log(4 * (sin((s(ii) - s(jj)) / 2))^2 / exp(1));
end

end