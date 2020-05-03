function res = x_vector(l,t,size)

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN x_vector(i)' newline 'l = ', num2str(l)]);
end

param = 3;

global q_m;

% Inner
if l == 1
    %res = radial(t, q_m, size) .* cos(t);
    res(1, 1:8) = (cos(t).^2 + 0.25 * sin(t).^2) .* cos(t);
% Outer
elseif l == 2
    res(1, 1:8) = 2 * cos(t);
else

%[x1, x2] = boundary(l, t);
%res = x1;

end