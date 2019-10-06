function res = r1(l, t)
% Calculates the length of the vector
% on a given boundary.
% As parameters takes:
%   l - boundary (1 or 2)
%   t - from [0, 2pi]

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN r(l, t)' newline 'l = ', num2str(l)]);
end

x = x_vector(l, t);
y = y_vector(l, t);

res = sqrt(x.^2 + y.^2);

end
