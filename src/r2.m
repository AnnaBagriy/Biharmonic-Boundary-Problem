function res = r2(l1, l2, t1, t2)
% Calculates the distance between two points
% on the given boundaries
% As parameters takes:
%   l1, l2 - boundaries (1 or 2)
%   t1, t2 are from [0, 2pi]

if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
    error(['WRONG INDEX IN r2(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
end

global x1;
global y1;
global x2;
global y2;

if l1 == 1
    if l2 == 1
        if t1 == t2
           res = 0; 
        else
            res = sqrt((x1(t1) - x1(t2))^2 + (y1(t1) - y1(t2))^2);
        end
    else
        res = sqrt((x1(t1) - x2(t2))^2 + (y1(t1) - y2(t2))^2);
    end
else
    if l2 == 1
        res = sqrt((x2(t1) - x1(t2))^2 + (y2(t1) - y1(t2))^2);
    else
        if t1 == t2
           res = 0; 
        else
           res = sqrt((x2(t1) - x2(t2))^2 + (y2(t1) - y2(t2))^2);
        end
    end
end
%res
end