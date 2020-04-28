classdef Calculator < handle
    
    properties (SetAccess = public)
       
       x1
       y1
       
       x2
       y2
       
       n1_on_1
       n2_on_1
       
       n1_on_2
       n2_on_2
       
    end
    
    methods (Access = public)
       
        function obj = Calculator(x1, x2, y1, y2, n1_on_1, n2_on_1, n1_on_2, n2_on_2)
            obj.x1 = x1;
            obj.x2 = x2;
            obj.y1 = y1;
            obj.y2 = y2;
            
            obj.n1_on_1 = n1_on_1;
            obj.n2_on_1 = n2_on_1;
            obj.n1_on_2 = n1_on_2;
            obj.n2_on_2 = n2_on_2;
        end
        
        function res = r2(obj, l1, l2, t1, t2)
        % Calculates the distance between two points
        % on the given boundaries
        % As parameters takes:
        %   l1, l2 - boundaries (1 or 2)
        %   t1, t2 are from [0, 2pi]

        if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
            error(['WRONG INDEX IN r2(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        end

        if l1 == 1
            if l2 == 1
                if t1 == t2
                   res = 0; 
                else
                    res = sqrt((obj.x1(t1) - obj.x1(t2)).^2 + (obj.y1(t1) - obj.y1(t2)).^2);
                end
            else
                res = sqrt((obj.x1(t1) - obj.x2(t2)).^2 + (obj.y1(t1) - obj.y2(t2)).^2);
            end
        else
            if l2 == 1
                res = sqrt((obj.x2(t1) - obj.x1(t2)).^2 + (obj.y2(t1) - obj.y1(t2)).^2);
            else
                if t1 == t2
                   res = 0; 
                else
                   res = sqrt((obj.x2(t1) - obj.x2(t2)).^2 + (obj.y2(t1) - obj.y2(t2)).^2);
                end
            end
        end
        %res
        end
        
        function res = H1(obj, l1, l2, ii, jj)
        % First type of H function
        % As parameters takes:
        %   l1 - one of the boundaries (1 or 2)
        %   l2 - one of the boundaries (1 or 2)
        %   l1 != l2
        %   t1, t2 - are from [0, 2pi]

        if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
            error(['WRONG INDEX IN H1(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        elseif l1 == l2
            error(['l1 CANNOT BE EQUAL TO l2 IN H1(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        end

        res = obj.r2(l1, l2, ii, jj)^2 * log(obj.r2(l1, l2, ii, jj)) / 4;

        end
        
    end
    
end