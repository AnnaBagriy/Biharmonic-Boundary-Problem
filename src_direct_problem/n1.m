function res = n1(l, tt)
% First component of the norm vector
% As parameters takes:
%   l - boundary (1 or 2)
%   t - is from [0, 2pi]

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN n1(l, t)' newline 'l = ', num2str(l)]);
end
syms t;

y_derivative = matlabFunction(diff(y_vector(l,t),t));
x_derivative = matlabFunction(diff(x_vector(l,t),t));

res = y_derivative(tt) ./ sqrt(x_derivative(tt).^2 + y_derivative(tt).^2);

end
