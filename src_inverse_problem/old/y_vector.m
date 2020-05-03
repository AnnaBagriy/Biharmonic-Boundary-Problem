function res = y_vector(l,t,size)

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN y_vector(i)' newline 'l = ', num2str(l)]);
end

param = 3;

global q_m;

% Inner
if l == 1
    %res = radial(t, q_m, size) .* sin(t);
    res(1, 1:8) = (cos(t).^2 + 0.25 * sin(t).^2) .* sin(t);
% Outer
elseif l == 2
    res(1, 1:8) = 1.5 * sin(t);
else
    
% [x1, x2] = boundary(l, t);
% res = x2;

end
