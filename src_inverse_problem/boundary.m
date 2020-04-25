function [x_1, x_2] = boundary(l, s)
param = 3;
% TODO: introduce radius not parameter

global r_1;
global r_2;

% Inner
if l == 1
    x_1 = r_1;
    x_2 = r_2;
% Outer
elseif l == 2
    x_1 = param * cos(s);
    x_2 = param * sin(s);
else
    error('WRONG INDEX IN boundary(l, s)');
end

end