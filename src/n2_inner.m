function res = n2_inner(r, tt)
syms t;

y_derivative = matlabFunction(diff(y_vector_inner(r,t),t));
x_derivative = matlabFunction(diff(x_vector_inner(r,t),t));

res = - x_derivative(tt) ./ sqrt(x_derivative(tt).^2 + y_derivative(tt).^2);

end
