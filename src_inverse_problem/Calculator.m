classdef Calculator < handle
    
    properties (SetAccess = public)
       
       s
        
       x1
       y1
       
       x2
       y2
       
       n1_on_1
       n2_on_1
       
       n1_on_2
       n2_on_2
       
       x_derivative_1
       y_derivative_1
      
       x_derivative_2
       y_derivative_2
       
      x_second_derivative_1
      y_second_derivative_1
      
      x_second_derivative_2
      y_second_derivative_2
      
    end
     
    properties (SetAccess = private)
       
       m
       
    end
    
    methods (Access = public)
       
        function obj = Calculator(s)
            obj.m = size(s, 2) / 2;

            obj.s = s;
        end
        
        %---------------------%
        
        function SetBoundary1(obj, x1, y1, q_m)
            obj.x1 = x1 + q_m;
            obj.y1 = y1;
        end
        
        function SetBoundary2(obj, x2, y2)
            obj.x2 = x2;
            obj.y2 = y2;
        end
        
        function SetNormal1(obj, n1, n2)
            obj.n1_on_1 = n1;
            obj.n2_on_1 = n2;
        end
        
        function SetNormal2(obj, n1, n2)
            obj.n1_on_2 = n1;
            obj.n2_on_2 = n2;
        end
        
        function SetDerivative1(obj, x_derivative_1, y_derivative_1)
            obj.x_derivative_1 = x_derivative_1;
            obj.y_derivative_1 = y_derivative_1;
        end
        
        function SetDerivative2(obj, x_derivative_2, y_derivative_2)
            obj.x_derivative_2 = x_derivative_2;
            obj.y_derivative_2 = y_derivative_2;
        end
        
        function SetSecondDerivative1(obj, x_2_derivative_1, y_2_derivative_1)
            obj.x_second_derivative_1 = x_2_derivative_1;
            obj.y_second_derivative_1 = y_2_derivative_1;
        end
        
        function SetSecondDerivative2(obj, x_2_derivative_2, y_2_derivative_2)
            obj.x_second_derivative_2 = x_2_derivative_2;
            obj.y_second_derivative_2 = y_2_derivative_2;
        end
        
        %---------------------%
        
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
        
        function res = H1_1(obj, l, ii, jj)
        % First type of H(1) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H1_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        res = r2(l, l, ii, jj).^2 ./ 8;

        end
        
        function res = H1_2(obj, l, ii, jj)
        % First type of H(2) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H1_2(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if ii == jj
            res = 0;
        else
            res = r2(l, l, ii, jj).^2 .* log((exp(1) .* r2(l, l, ii, jj)^2) ./ (4 .* (sin((obj.s(ii) - obj.s(jj)) ./ 2)).^2)) ./ 8;
        end

        end
        
        function res = H2(obj, l1, l2, ii, jj)
        % Second type of H function
        % As parameters takes:
        %   l1 - one of the boundaries (1 or 2)
        %   l2 - one of the boundaries (1 or 2)
        %   l1 != l2
        %   t1, t2 - are from [0, 2pi]

        if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
            error(['WRONG INDEX IN H2(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        elseif l1 == l2
            error(['l1 CANNOT BE EQUAL TO l2 IN H2(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        end

        if l2 == 1
            res = - (obj.n1_on_1(jj) .* (obj.x2(ii) - obj.x1(jj)) .* (1 + 2 .* log(r2(l1, l2, ii, jj))) +...
                     obj.n2_on_1(jj) .* (obj.y2(ii) - obj.y1(jj)) .* (1 + 2 .* log(r2(l1, l2, ii, jj)))) ./ 4;
        else
            res = - (obj.n1_on_2(jj) .* (obj.x1(ii) - obj.x2(jj)) .* (1 + 2 .* log(r2(l1, l2, ii, jj))) +...
                     obj.n2_on_2(jj) .* (obj.y1(ii) - obj.y2(jj)) .* (1 + 2 .* log(r2(l1, l2, ii, jj)))) ./ 4;
        end

        end
        
        function res = H2_1(obj, l, ii, jj)
        % Second type of H(1) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector 
        %   jj - index for s-vector 
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H2_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if ii == jj
            res = 0;
        else
            if l == 1
                res = -(obj.n1_on_1(jj) * (obj.x1(ii) - obj.x1(jj)) + obj.n2_on_1(jj) * (obj.y1(ii) - obj.y1(jj))) / 4;
            else
                res = -(obj.n1_on_2(jj) * (obj.x2(ii) - obj.x2(jj)) + obj.n2_on_2(jj) * (obj.y2(ii) - obj.y2(jj))) / 4;
            end
        end

        end
        
        function res = H2_2(obj, l, ii, jj)
        % Second type of H(2) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H2_2(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if ii == jj 
            res = 0;
        else
            if l == 1
                res = (obj.n1_on_1(jj) * (obj.x1(ii) - obj.x1(jj)) * (log((4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2) / (exp(1) * r2(l, l, ii, jj)^2)) - 1) +...
                       obj.n2_on_1(jj) * (obj.y1(ii) - obj.y1(jj)) * (log((4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2) / (exp(1) * r2(l, l, ii, jj)^2)) - 1)) / 4;
            else
                res = (obj.n1_on_2(jj) * (obj.x2(ii) - obj.x2(jj)) * (log((4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2) / (exp(1) * r2(l, l, ii, jj)^2)) - 1) +...
                       obj.n2_on_2(jj) * (obj.y2(ii) - obj.y2(jj)) * (log((4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2) / (exp(1) * r2(l, l, ii, jj)^2)) - 1)) / 4;
            end    
        end

        end
        
        function res = H3(obj, l1, l2, ii, jj)
        % Third type of H function
        % As parameters takes:
        %   l1 - one of the boundaries (1 or 2)
        %   l2 - one of the boundaries (1 or 2)
        %   l1 != l2
        %   t1, t2 - are from [0, 2pi]

        if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
            error(['WRONG INDEX IN H3(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        elseif l1 == l2
            error(['l1 CANNOT BE EQUAL TO l2 IN H3(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        end

        if l1 == 1
            res = (obj.n1_on_1(ii) * (obj.x1(ii) - obj.x2(jj)) * (1 + 2 * log(r2(l1, l2, ii, jj))) +...
                   obj.n2_on_1(ii) * (obj.y1(ii) - obj.y2(jj)) * (1 + 2 * log(r2(l1, l2, (ii), jj)))) / 4;
        else
            res = (obj.n1_on_2(ii) * (obj.x2(ii) - obj.x1(jj)) * (1 + 2 * log(r2(l1, l2, (ii), (jj)))) +...
                   obj.n2_on_2(ii) * (obj.y2(ii) - obj.y1(jj)) * (1 + 2 * log(r2(l1, l2, (ii), (jj))))) / 4;
        end

        end
        
        function res = H3_1(obj, l, ii, jj)
        % Third type of H(1) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H3_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if ii == jj
            res = 0;
        else
            if l == 1
                res = (obj.n1_on_1(ii) * (obj.x1(ii) - obj.x1(jj)) + obj.n2_on_1(ii) * (obj.y1(ii) - obj.y1(jj))) / 4;
            else
                res = (obj.n1_on_2(ii) * (obj.x2(ii) - obj.x2(jj)) + obj.n2_on_2(ii) * (obj.y2(ii) - obj.y2(jj))) / 4;
            end
        end

        end
        
        function res = H3_2(obj, l, ii, jj)
        % Third type of H(2) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H3_2(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if ii == jj
            res = 0;
        else
            if l == 1
                    res = (obj.n1_on_1(ii) * (obj.x1(ii) - obj.x1(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2)) + 1)... 
                         + obj.n2_on_1(ii) * (obj.y1(ii) - obj.y1(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2)) + 1)) / 4;
            else
                    res = (obj.n1_on_2(ii) * (obj.x2(ii) - obj.x2(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2)) + 1)... 
                         + obj.n2_on_2(ii) * (obj.y2(ii) - obj.y2(jj)) * (log((exp(1) * r2(l, l, (ii), (jj))^2) / (4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2)) + 1)) / 4;
            end
        end

        end
        
        function res = H4(obj, l1, l2, ii, jj)
        % Forth type of H function
        % As parameters takes:
        %   l1 - one of the boundaries (1 or 2)
        %   l2 - one of the boundaries (1 or 2)
        %   l1 != l2
        %   t1, t2 - are from [0, 2pi]

        %if l1 ~= 1 && l1 ~= 2 && l2 ~= 1 && l2 ~= 2
        %    error(['WRONG INDEX IN H4(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        %elseif l1 == l2
        %    error(['l1 CANNOT BE EQUAL TO l2 IN H4(l1, l2, t1, t2)' newline 'l1 = ', num2str(l1) newline 'l2 = ', num2str(l2)]);
        %end

        if l1 == 1 && l2 == 1
            n_1c = obj.n1_on_1(ii);
            n_2c = obj.n2_on_1(ii);

            n_1i = obj.n1_on_1(jj);
            n_2i = obj.n2_on_1(jj);

            x_c = obj.x1(ii);
            y_c = obj.y1(ii);

            x_i = obj.x1(jj);
            y_i = obj.y1(jj);
        end
        if l1 == 1 && l2 == 2
            n_1c = obj.n1_on_1(ii);
            n_2c = obj.n2_on_1(ii);

            n_1i = obj.n1_on_2(jj);
            n_2i = obj.n2_on_2(jj);

            x_c = obj.x1(ii);
            y_c = obj.y1(ii);

            x_i = obj.x2(jj);
            y_i = obj.y2(jj);
        end
        if l1 == 2 && l2 == 1
            n_1c = obj.n1_on_2(ii);
            n_2c = obj.n2_on_2(ii);

            n_1i = obj.n1_on_1(jj);
            n_2i = obj.n2_on_1(jj);

            x_c = obj.x2(ii);
            y_c = obj.y2(ii);

            x_i = obj.x1(jj);
            y_i = obj.y1(jj);
        end
        if l1 == 2 && l2 == 2
            n_1c = obj.n1_on_2(ii);
            n_2c = obj.n2_on_2(ii);

            n_1i = obj.n1_on_2(jj);
            n_2i = obj.n2_on_2(jj);

            x_c = obj.x2(ii);
            y_c = obj.y2(ii);

            x_i = obj.x2(jj);
            y_i = obj.y2(jj);
        end

        res = - (n_1c * (x_c - x_i) + n_2c * (y_c - y_i)) * (n_1i * (x_c - x_i) + n_2i * (y_c - y_i)) /...
                 2 / r2(l1, l2, ii, jj)^2 - (n_1c * n_1i + n_2c * n_2i) * (1 + 2 * log(r2(l1, l2, ii, jj))) / 4;
        end
        
        function res = H4_1(obj, l, ii, jj)
        % Fourth type of H(1) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if l == 1
            res = - (obj.n1_on_1(ii) * obj.n1_on_1(jj) + obj.n2_on_1(ii) * obj.n2_on_1(jj)) / 4;
        else
            res = - (obj.n1_on_2(ii) * obj.n1_on_2(jj) + obj.n2_on_2(ii) * obj.n2_on_2(jj)) / 4;  
        end

        end
        
        function res = H4_2(obj, l, ii, jj)
        % Fourth type of H(2) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if ii == jj
            if l == 1
                x_der = obj.x_derivative_1(ii);
                y_der = obj.y_derivative_1(ii);
            else
                x_der = obj.x_derivative_2(ii);
                y_der = obj.y_derivative_2(ii);
            end

            res = - (1 + log(exp(1) * ((y_der).^2 + (x_der).^2))) / 4;
        else
            res = obj.H4(l,l, ii, jj) - obj.H4_1(l, ii, jj) .* log(4 * (sin((obj.s(ii) - obj.s(jj)) / 2)).^2 / exp(1));
        end

        end
        
        function res = H5(obj, l1, l2, ii, jj, nu)
        % Fifth type of H function
        % As parameters takes:
        %   l1 - one of the boundaries (1 or 2)
        %   l2 - one of the boundaries (1 or 2)
        %   l1 != l2
        %   t1, t2 - are from [0, 2pi]

        if l1 == 1 && l2 == 1
            n_1c = obj.n1_on_1(ii);
            n_2c = obj.n2_on_1(ii);

            x_c = obj.x1(ii);
            y_c = obj.y1(ii);

            x_i = obj.x1(jj);
            y_i = obj.y1(jj);
        end
        if l1 == 1 && l2 == 2
            n_1c = obj.n1_on_1(ii);
            n_2c = obj.n2_on_1(ii);

            x_c = obj.x1(ii);
            y_c = obj.y1(ii);

            x_i = obj.x2(jj);
            y_i = obj.y2(jj);
        end
        if l1 == 2 && l2 == 1
            n_1c = obj.n1_on_2(ii);
            n_2c = obj.n2_on_2(ii);

            x_c = obj.x2(ii);
            y_c = obj.y2(ii);

            x_i = obj.x1(jj);
            y_i = obj.y1(jj);
        end
        if l1 == 2 && l2 == 2
            n_1c = obj.n1_on_2(ii);
            n_2c = obj.n2_on_2(ii);

            x_c = obj.x2(ii);
            y_c = obj.y2(ii);

            x_i = obj.x2(jj);
            y_i = obj.y2(jj);
        end

        res = (1 + 3 * nu) / 8 * pi + ((1 - nu) * (n_1c * (x_c - x_i) +...
            n_2c * (y_c - y_i))^2) / ( 4 * pi * r2(l1, l2, ii, jj)^2) +...
            ((1 + nu) * log(r2(l1, l2, ii, jj)^2)) / (8 * pi);

        end
        
        function res = H5_1(obj, l, ii, jj, nu)
        % Fifth type of H(1) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        res = (1 + nu) / 4;

        end
        
        function res = H5_2(obj, l, ii, jj, nu)
        % Fifth type of H(2) function
        % As parameters takes:
        %   l - boundaries (1 or 2)
        %   l1 = l2 = l
        %   ii - index for s-vector
        %   jj - index for s-vector
        %   s - vector

        if l ~= 1 && l ~= 2
            error(['WRONG INDEX IN H4_1(l, ii, jj, s)' newline 'l = ', num2str(l)]);
        end

        if l == 1
             if ii == jj
                 res = (1 + 3 * nu) / 4 + (1 + nu) * log(exp(1) * (obj.x_derivative_1(ii)^2 + obj.y_derivative_1(ii)^2)^2) / 4;
             else
                 res = obj.H5(l, l, ii, jj, nu) - (1 + nu) * log(4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2  / exp(1)) / 4;
             end
        else
             if ii == jj
                 res = (1 + 3 * nu) / 4 + (1 + nu) * log(exp(1) * (obj.x_derivative_2(ii)^2 + obj.y_derivative_2(ii)^2)^2) / 4;
             else
                 res = obj.H5(l, l, ii, jj, nu) - (1 + nu) * log(4 * (sin((obj.s(ii) - obj.s(jj)) / 2))^2  / exp(1)) / 4;
             end
        end

        end
        
        function res = H6(obj, l1, l2, ii, jj, nu)
        % Sixth type of H function
        % As parameters takes:
        %   l1 - one of the boundaries (1 or 2)
        %   l2 - one of the boundaries (1 or 2)
        %   t1, t2 - are from [0, 2pi]

        if l1 == 1 && l2 == 1
            n_1c = obj.n1_on_1(ii);
            n_2c = obj.n2_on_1(ii);

            n_1i = obj.n1_on_1(jj);
            n_2i = obj.n2_on_1(jj);

            x_c = obj.x1(ii);
            y_c = obj.y1(ii);

            x_i = obj.x1(jj);
            y_i = obj.y1(jj);

            x_1der_1 = obj.x_derivative_1(ii);
            x_1der_2 = obj.y_derivative_1(ii);

            x_2der_1 = obj.x_second_derivative_1(ii);
            x_2der_2 = obj.y_second_derivative_1(ii);
        end
        if l1 == 1 && l2 == 2
            n_1c = obj.n1_on_1(ii);
            n_2c = obj.n2_on_1(ii);

            n_1i = obj.n1_on_2(jj);
            n_2i = obj.n2_on_2(jj);

            x_c = obj.x1(ii);
            y_c = obj.y1(ii);

            x_i = obj.x2(jj);
            y_i = obj.y2(jj);

            x_1der_1 = obj.x_derivative_1(ii);
            x_1der_2 = obj.y_derivative_1(ii);

            x_2der_1 = obj.x_second_derivative_1(ii);
            x_2der_2 = obj.y_second_derivative_1(ii);
        end
        if l1 == 2 && l2 == 1
            n_1c = obj.n1_on_2(ii);
            n_2c = obj.n2_on_2(ii);

            n_1i = obj.n1_on_1(jj);
            n_2i = obj.n2_on_1(jj);

            x_c = obj.x2(ii);
            y_c = obj.y2(ii);

            x_i = obj.x1(jj);
            y_i = obj.y1(jj);

            x_1der_1 = obj.x_derivative_2(ii);
            x_1der_2 = obj.y_derivative_2(ii);

            x_2der_1 = obj.x_second_derivative_2(ii);
            x_2der_2 = obj.y_second_derivative_1(ii);
        end
        if l1 == 2 && l2 == 2
            n_1c = obj.n1_on_2(ii);
            n_2c = obj.n2_on_2(ii);

            n_1i = obj.n1_on_2(jj);
            n_2i = obj.n2_on_2(jj);

            x_c = obj.x2(ii);
            y_c = obj.y2(ii);

            x_i = obj.x2(jj);
            y_i = obj.y2(jj);

            x_1der_1 = obj.x_derivative_2(ii);
            x_1der_2 = obj.y_derivative_2(ii);

            x_2der_1 = obj.x_second_derivative_2(ii);
            x_2der_2 = obj.y_second_derivative_2(ii);
        end

        if ii == jj
            res = (1 - 3 * nu) * (n_1c * x_2der_1 + n_2c * x_2der_2) /...
                (4 * (x_1der_1^2 + x_1der_2^2));
        else
            res = (1 - nu) * (((n_1c * (x_c - x_i) + n_2c *...
                (y_c - y_i))^2 * (n_1i * (x_c - x_i) +...
                n_2i * (y_c - y_i))) / r2(l1, l2, ii, jj)^4 -...
                ((n_1c * n_1i + n_2c * n_2i) * (n_1c *...
                (x_c - x_i) + n_2c * (y_c - y_i))) /...
                r2(l1, l2, ii, jj)^2) / (2 * pi) - ((1 + nu) * (n_1i *...
                (x_c - x_i) + n_2i * (y_c - y_i)) /...
                (4 * pi * r2(l1, l2, ii, jj)^2));
        end

        end
        
        function res = R(obj, k, m)
        % Calculates weight function
        % for Nystrom method

        jj = 1:m - 1;
        r = cos(jj .* pi .* k / m) ./ jj;

        res = - (0.5 + (-1)^k / (2 * m) + sum(r)) / m;

        end
        
        function res = N_r_1(obj, ii, jj)

        res = (cos(obj.s(jj)) * obj.n1_on_2(jj) + sin(obj.s(jj)) * obj.n2_on_2(jj)) * (log(r2(2, 1, ii, jj)) + 3);

        end

        function res = N_r_2(obj, ii, jj)

        res = (cos(obj.s(jj)) * (obj.x2(ii) - obj.x1(jj)) + sin(obj.s(jj)) * (obj.y2(ii) - obj.y1(jj))) * (log(r2(2, 1, ii, jj)) + 1);

        end
        
        function res = l(obj, ii, t, m)

        if ii <= m + 1
            res = cos((ii - 1) .* t);
        else
            res = sin(m - ii - 1) .* t;
        end

        end

    end
    
end