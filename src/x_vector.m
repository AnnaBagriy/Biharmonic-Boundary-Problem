function res = x_vector(l,t)

if l ~= 1 && l ~= 2
    error(['WRONG INDEX IN x_vector(i)' newline 'l = ', num2str(l)]);
end

[x1, x2] = boundary(l, t);
res = x1;

end