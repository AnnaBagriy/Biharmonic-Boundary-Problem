% Display area D
tn=0:0.01:2 * pi;

x1 = x_vector(1, tn);
x2 = x_vector(2, tn);
y1 = y_vector(1, tn);
y2 = y_vector(2, tn);

 patch([x1, x2],...
       [y1, y2],...
       [.7 .7 .7])

matrix_file = 'A_matrix.txt';
vector_file = 'y_vector.txt';

A = dlmread(matrix_file);
y = dlmread(vector_file);

m = (size(y, 1) - 3) / 8;

global s;

r = @(x1, y1, x2, y2) sqrt((x1 - x1).^2 + (y1 - y2).^2);
s11=-3/2;
s12=3/2;
s21=2;
s22=2;
s31=3;
s32=-3;
s41=-5/2;
s42=-5/2;
format long;
% Exact solution
% U = @(x,y) (10 * r(x,y,s11,s12)^2 * log(r(x,y,s11,s12)) + 4 * r(x,y,s21,s22)^2 * log(r(x,y,s21,s22)) + ...
%     3 * r(x,y,s21,s22)^2 * log(r(x,y,s21,s22)) + 2 * r(x,y,s41,s42)^3 * log(r(x,y,s41,s42))) / (8 * pi);

U = @(x,y) 1; %x-2.*y;

det_A = det(A);
x = A \ y;

a0 = x(8 * m + 1);
a1 = x(8 * m + 2);
a2 = x(8 * m + 3);

fi_func_1 = zeros(2 * m, 1);
fi_func_2 = zeros(2 * m, 1);
psi_func_1 = zeros(2 * m, 1);
psi_func_2 = zeros(2 * m, 1);

for ii = 1:2 * m 
    fi_func_1(ii) = x(ii);
    fi_func_2(ii) = x(ii + 2 * m);
    psi_func_1(ii) = x(ii + 4 * m);
    psi_func_2(ii) = x(ii + 6 * m);
end

% Find approximate solution U
U_approx = 0;

x = 1.5;
y = 0;
%  x=-0.75;
%  y=-1;

r_1 = sqrt((x_vector_inner(1, s) - x).^2 + (y_vector_inner(1, s) - y).^2);
r_2 = sqrt((x_vector_inner(2, s) - x).^2 + (y_vector_inner(2, s) - y).^2);

syms t;
y_derivative1 = matlabFunction(diff(y_vector_inner(1, t),t));
x_derivative1 = matlabFunction(diff(x_vector_inner(1, t),t));
y_derivative2 = matlabFunction(diff(y_vector_inner(2, t),t));
x_derivative2 = matlabFunction(diff(x_vector_inner(2, t),t));

H11 =  r_1.^2 .* log(r_1) ./ 4;
H21 =  -(n1_inner(1, s) .* (x - x_vector_inner(1, s)) .* (1 + 2 .* log(r_1)) +...
         n2_inner(1, s) .* (y - y_vector_inner(1, s)) .* (1 + 2 .* log(r_1))) ./ 4;
     
H12 =  r_2.^2 .* log(r_2) ./ 4;
H22 =  -(n1_inner(2, s) .* (x - x_vector_inner(2, s)) .* (1 + 2 .* log(r_2)) +...
         n2_inner(2, s) .* (y - y_vector_inner(2, s)) .* (1 + 2 .* log(r_2))) ./ 4;

for ii = 1:2 * m 
      U_approx = U_approx + fi_func_1(ii) *...
          H11(ii) + psi_func_1(ii) * H21(ii) +...
               fi_func_2(ii) * H12(ii) + psi_func_2(ii) * H22(ii);
end
U_approx = U_approx/2/m + a0 + a1 * x + a2 * y;
U_exact = U(x, y);

error = abs(U_approx - U_exact);

disp([newline char(9) 'RESULTS'])
disp(['m = ' num2str(m)])
disp(['U_exact('  num2str(x) ', ' num2str(y) ') = ' num2str(U_exact)])
disp(['U_approx(' num2str(x) ', ' num2str(y) ') = ' num2str(U_approx)])
disp(['error = ' num2str(error)])
