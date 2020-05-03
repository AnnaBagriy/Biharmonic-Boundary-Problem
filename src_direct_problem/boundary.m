function [x_1, x_2] = boundary(l, s)
param = 2;
% TODO: introduce radius not parameter

r = @(t) sqrt(cos(t).^2 + 0.25 * sin(t).^2);

% Inner
if l == 1
    x_1 = r(s) .* cos(s);
    x_2 = r(s) .* sin(s);
% Outer
elseif l == 2
    x_1 = param * cos(s);
    x_2 = param * sin(s);
else
    error('WRONG INDEX IN boundary(l, s)');
end

end