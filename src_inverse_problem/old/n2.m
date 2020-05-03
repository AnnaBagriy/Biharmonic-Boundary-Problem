function res = n2(l,tt)
% Second component of the norm vector
% As parameters takes:
%   l - boundary (1 or 2)
%   tt - is from [0, 2pi]

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN n2(l, t1, t2)' newline 'l = ', num2str(l)]);
end

syms t;

y_derivative = matlabFunction(diff(y_vector(l,t, size(t, 2)),t));
x_derivative = matlabFunction(diff(x_vector(l,t, size(t, 2)),t));

x_derivative_calc = zeros(1, size(tt, 2));
y_derivative_calc = zeros(1, size(tt, 2));

for ii = 1:size(tt, 2)
    x = x_derivative(tt(ii));
    y = y_derivative(tt(ii));

    x_derivative_calc(ii) = x(ii);
    y_derivative_calc(ii) = y(ii);
end

res = - x_derivative_calc ./ sqrt(x_derivative_calc.^2 + y_derivative_calc.^2);

end