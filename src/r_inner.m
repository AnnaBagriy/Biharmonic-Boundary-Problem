function res = r_inner(r1, r2, t1, t2)
% Calculates the distance between two points
% on the given boundaries
% As parameters takes:
%   l1, l2 - boundaries (1 or 2)
%   t1, t2 are from [0, 2pi]

res = sqrt((x_vector_inner(r1, t1) - x_vector_inner(r2, t2))^2 + (y_vector_inner(r1, t1) - y_vector_inner(r2, t2))^2);

end