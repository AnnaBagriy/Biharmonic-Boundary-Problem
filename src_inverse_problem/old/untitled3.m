
classdef MainSystem
   properties
      Value {mustBeNumeric}
      m
      h = pi / m
      s = (1:2 * m - 1) .* h
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value] * n;
      end
   end
end
m = 4;
h = pi / m;
jj = 1:2 * m;
s = (jj - 1) * h;

x_on_1 = cos(s);
y_on_1 = sin(s);
x_on_2 = 2 * cos(s);
y_on_2 = 2 * sin(s);

n1_on_1 = n1(1, s);
n2_on_1 = n2(1, s);
n1_on_2 = n1(2, s);
n2_on_2 = n2(2, s);