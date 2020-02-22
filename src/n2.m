function res = n2(l,tt)
% Second component of the norm vector
% As parameters takes:
%   l - boundary (1 or 2)
%   tt - is from [0, 2pi]

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN n2(l, t1, t2)' newline 'l = ', num2str(l)]);
end

syms t;

y_derivative = matlabFunction(diff(y_vector(l,t),t));
x_derivative = matlabFunction(diff(x_vector(l,t),t));

res = - x_derivative(tt) ./ sqrt(x_derivative(tt).^2 + y_derivative(tt).^2);

end