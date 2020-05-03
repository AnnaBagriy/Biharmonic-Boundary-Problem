global s;
global n_11; % 1 compoment on l=1
global n_21; % 2 compoment on l=1
global n_12; % 1 compoment on l=2
global n_22; % 2 compoment on l=2
global x1;
global y1;
global x2;
global y2;
global m;
format long;
m = 16;
h = pi / m;
jj = 1:2 * m;
s = (jj - 1) * h;

r_1 = @(x,y) sqrt((x_vector(1, s) - x).^2 + (y_vector(1, s) - y).^2);
r_2 = @(x,y) sqrt((x_vector(2, s) - x).^2 + (y_vector(2, s) - y).^2);

% Utility coeffs
A0 = 1;
A1 = 1;
A2 = 1;

% Main matrix
A = zeros(8 * m + 3, 8 * m + 3);
y = zeros (8 * m + 3, 1);

% Right side of the system
y(1) = A0;
y(2) = A1;
y(3) = A2;

x1 = x_vector(1, s);
y1 = y_vector(1, s);
x2 = x_vector(2, s);
y2 = y_vector(2, s);

n_11 = n1(1, s);
n_21 = n2(1, s);
n_12 = n1(2, s);
n_22 = n2(2, s);

% Exact solution
U = @(x,y) (10 * r_1(s11,s12)^2 * log(r_1(s11,s12)) + 4 * r_1(s21,s22)^2 * log(r_1(s21,s22)) + ...
    3 * r_1(s31,s32)^2 * log(r_1(s31,s32)) + 2 * r_1(s41,s42)^3 * log(r_1(s41,s42))) / (8 * pi);

% Define functions on boundaries
% f1 = @(x, y) (10 .* r_1(s11,s12).^2 .* log(r_1(s11,s12)) + 4 * r_1(s21,s22).^2 .* log(r_1(s21,s22)) + ...
%     3 * r_1(s31,s32).^2 .* log(r_1(s31,s32)) + 2 * r_1(s41,s42).^3 .* log(r_1(s41,s42))) / (8 * pi);
% f2 = @(x, y) (10 * r_2(s11,s12).^2 .* log(r_2(s11,s12)) + 4 .* r_2(s21,s22).^2 .* log(r_2(s21,s22)) + ...
%     3 * r_2(s31,s32).^2 .* log(r_2(s31,s32)) + 2 .* r_2(s41,s42).^3 .* log(r_2(s41,s42))) / (8 * pi);

syms t;
y_derivative1 = matlabFunction(diff(y_vector(1, t),t));
x_derivative1 = matlabFunction(diff(x_vector(1, t),t));
y_derivative2 = matlabFunction(diff(y_vector(2, t),t));
x_derivative2 = matlabFunction(diff(x_vector(2, t),t));

f1 = @(x, y) x-2.*y;
f2 = @(x, y) x-2.*y;
g1 = @(x, y) n_11;
g2 = @(x, y) n_12;

f1_k = f1(x1, y1);
f2_k = f2(x1, y1);
g1_k = g1(x2, y2);
g2_k = g2(x2, y2);

%k = 4;
%k_loop = k;
for ii = 1:2 * m
    y(ii + 3) =1; %f1(x1(ii), y1(ii));
    y(2 * m + 3 + ii) = 1;%f2(x2(ii), y2(ii));
    y(4 * m + 3 + ii) = 0; %n_11(ii) -2* n_21(ii);
    y(6 * m + 3 + ii) =  0;%n_12(ii) -2* n_22(ii);
end

for ii = 1:2 * m   
    A(1, ii) = h;
    A(1, 2 * m + ii) = h;
    
    A(2, ii) = h * x1(ii);
    A(2, 2 * m + ii) = h * x2(ii);
    A(2, 4 * m + ii) = h * n_11(ii);
    A(2, 6 * m + ii) = h * n_12(ii);
    
    A(3, ii) = h * y1(ii);
    A(3, 2 * m + ii) = h * y2(ii);
    A(3, 4 * m + ii) = h * n_21(ii);
    A(3, 6 * m + ii) = h * n_22(ii);
    
    A(3 + ii, 8 * m + 1) = 1;
    A(3 + ii, 8 * m + 2) = x1(ii);
    A(3 + ii, 8 * m + 3) = y1(ii);
    
    A(3 + 2 * m + ii, 8 * m + 1) = 1;
    A(3 + 2 * m + ii, 8 * m + 2) = x2(ii);
    A(3 + 2 * m + ii, 8 * m + 3) = y2(ii);
    
    A(3 + 4 * m + ii, 8 * m + 1) = 0;
    A(3 + 4 * m + ii, 8 * m + 2) = n_11(ii);
    A(3 + 4 * m + ii, 8 * m + 3) = n_21(ii);
    
    A(3 + 6 * m + ii, 8 * m + 1) = 0;
    A(3 + 6 * m + ii, 8 * m + 2) = n_12(ii);
    A(3 + 6 * m + ii, 8 * m + 3) = n_22(ii);   

    for jj = 1:2 * m
        A(3 + ii, jj) = H1_1(1, ii, jj) * R(abs(jj - ii) , m) + H1_2(1, ii, jj) / (2 * m);
        A(3 + ii, jj + 2 * m) = H1(1, 2, ii, jj) / (2 * m);
        A(3 + ii, jj + 4 * m) = H2_1(1, ii, jj) * R(abs(jj - ii) , m) + H2_2(1, ii, jj) / (2 * m);
        A(3 + ii, jj + 6 * m) = H2(1, 2, ii, jj) / (2 * m);
        
        A(3 + ii + 2 * m, jj) = H1(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 2 * m, jj + 2 * m) = H1_1(2, ii, jj) * R(abs(jj - ii), m) + H1_2(2, ii, jj) / (2 * m);
        A(3 + ii + 2 * m, jj + 4 * m) = H2(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 2 * m, jj + 6 * m) = H2_1(2, ii, jj) * R(abs(jj - ii) , m) + H2_2(2, ii, jj) / (2 * m);
   
        A(3 + ii + 4 * m, jj) = H3_1(1, ii, jj) * R(abs(jj - ii) , m) + H3_2(1, ii, jj) / (2 * m);
        A(3 + ii + 4 * m, jj + 2 * m) = H3(1, 2, ii, jj) / (2 * m);
        A(3 + ii + 4 * m, jj + 4 * m) = H4_1(1, ii, jj) * R(abs(jj - ii), m) + H4_2(1, ii, jj) / (2 * m);
        A(3 + ii + 4 * m, jj + 6 * m) = H4(1, 2, ii, jj) / (2 * m);
        
        A(3 + ii + 6 * m, jj) = H3(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 6 * m, jj + 2 * m) = H3_1(2, ii, jj) * R(abs(jj - ii) , m) + H3_2(2, ii, jj) / (2 * m);
        A(3 + ii + 6 * m, jj + 4 * m) = H4(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 6 * m, jj + 6 * m) = H4_1(2, ii, jj) * R(abs(jj - ii), m) + H4_2(2, ii, jj) / (2 * m);
    end
end
    
A
y

matrix_file = 'A_matrix.txt';
vector_file = 'y_vector.txt';

mid = fopen(matrix_file, 'wt');
vid = fopen(vector_file, 'wt');

for ii = 1:size(A, 1)
    fprintf(mid, '%g\t', A(ii, :));
    fprintf(mid, '\n');
end

fprintf(vid, '%d\n', y);

fclose(mid);
fclose(vid);

%disp('the end of matrix initialization');
main