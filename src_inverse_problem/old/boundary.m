function [x_1, x_2] = boundary(l, s)
param = 3;
% TODO: introduce radius not parameter

global q_m;

% Inner
if l == 1
    x_1 = radial(s, q_m) .* cos(s);
    x_2 = radial(s, q_m) .* sin(s);
% Outer
elseif l == 2
    x_1 = param * cos(s);
    x_2 = param * sin(s);
else
    error('WRONG INDEX IN boundary(l, s)');
end

end