format long;

%---------------------------%

% Globals

global s;

global n_11; % 1 component on l = 1
global n_21; % 2 component on l = 1
global n_12; % 1 component on l = 2
global n_22; % 2 component on l = 2

global x1;
global y1;
global x2;
global y2;

global x_derivative1;
global x_derivative2;
global y_derivative1;
global y_derivative2;

global x_2derivative1;
global x_2derivative2;
global y_2derivative1;
global y_2derivative2;

global m;

global nu;

%---------------------------%

m = 8;
h = pi / m;
jj = 1:2 * m;
s = (jj - 1) * h;

nu = 0.5;

x1 = x_vector(1, s);
y1 = y_vector(1, s);
x2 = x_vector(2, s);
y2 = y_vector(2, s);

n_11 = n1(1, s);
n_21 = n2(1, s);
n_12 = n1(2, s);
n_22 = n2(2, s);

% Utility coeffs
A0 = 1;
A1 = 1;
A2 = 1;

%---------------------------%

% Derivatives

syms t;

y_derivative1_sym = matlabFunction(diff(y_vector(1, t), t));
y_derivative2_sym = matlabFunction(diff(y_vector(2, t), t));

y_2derivative1_sym = matlabFunction(diff(diff(y_vector(1, t), t), t));
y_2derivative2_sym = matlabFunction(diff(diff(y_vector(2, t), t), t));

x_derivative1_sym = matlabFunction(diff(x_vector(1, t), t));
x_derivative2_sym = matlabFunction(diff(x_vector(2, t), t));

x_2derivative1_sym = matlabFunction(diff(diff(x_vector(1, t), t), t));
x_2derivative2_sym = matlabFunction(diff(diff(x_vector(2, t), t), t));

y_derivative1 = y_derivative1_sym(s);
y_derivative2 = y_derivative2_sym(s);

x_derivative1 = x_derivative1_sym(s);
x_derivative2 = x_derivative2_sym(s);

x_2derivative1 = x_2derivative1_sym(s);
x_2derivative2 = x_2derivative2_sym(s);

y_2derivative1 = y_2derivative1_sym(s);
y_2derivative2 = y_2derivative2_sym(s);

%---------------------------%

% Right vector functions
g2 = @(x, y) x + y;
q2 = @(x, y) x - y;

g2_k = g2(x2, y2);
q2_k = q2(x2, y2);

% Right vector initialization
y = zeros (8 * m + 3, 1);

y(1) = A0;
y(2) = A1;
y(3) = A2;

for ii = 1:2 * m
    y(ii + 3) = 0;
    y(2 * m + 3 + ii) = 0;
    y(4 * m + 3 + ii) = g2_k(ii);
    y(6 * m + 3 + ii) =  q2_k(ii);
end

% Matrix initialization
A = zeros(8 * m + 3, 8 * m + 3);

for ii = 1:2 * m   
    
    %---------------------------%
    
    % 1st utility equation
    A(1, ii) = h;
    A(1, 2 * m + ii) = h;
    
    % 2nd utility equation
    A(2, ii) = h * x1(ii);
    A(2, 2 * m + ii) = h * x2(ii);
    A(2, 4 * m + ii) = h * n_11(ii);
    A(2, 6 * m + ii) = h * n_12(ii);
    
    % 3d utility equation
    A(3, ii) = h * y1(ii);
    A(3, 2 * m + ii) = h * y2(ii);
    A(3, 4 * m + ii) = h * n_21(ii);
    A(3, 6 * m + ii) = h * n_22(ii);
    
    %---------------------------%
    
    % 1st equation linear function
    A(3 + ii, 8 * m + 1) = 1;
    A(3 + ii, 8 * m + 2) = x1(ii);
    A(3 + ii, 8 * m + 3) = y1(ii);
    
    % 2nd equation linear function
    A(3 + 2 * m + ii, 8 * m + 1) = 0;
    A(3 + 2 * m + ii, 8 * m + 2) = n_11(ii);
    A(3 + 2 * m + ii, 8 * m + 3) = n_21(ii);
    
    % 3d equation linear function
    A(3 + 4 * m + ii, 8 * m + 1) = 0;
    A(3 + 4 * m + ii, 8 * m + 2) = n_12(ii);
    A(3 + 4 * m + ii, 8 * m + 3) = n_22(ii);
    
    % 4th equation linear function
    A(3 + 6 * m + ii, 8 * m + 1) = 0;
    A(3 + 6 * m + ii, 8 * m + 2) = 0;
    A(3 + 6 * m + ii, 8 * m + 3) = 0;   

    %---------------------------%
    
    for jj = 1:2 * m
        % 1st equation kernels
        A(3 + ii, jj) = H1_1(1, ii, jj) * R(abs(jj - ii) , m) + H1_2(1, ii, jj) / (2 * m);
        A(3 + ii, jj + 2 * m) = H1(1, 2, ii, jj) / (2 * m);
        A(3 + ii, jj + 4 * m) = H2_1(1, ii, jj) * R(abs(jj - ii) , m) + H2_2(1, ii, jj) / (2 * m);
        A(3 + ii, jj + 6 * m) = H2(1, 2, ii, jj) / (2 * m);
        
        % 2nd equation kernels
        A(3 + ii + 2 * m, jj) = H3_1(1, ii, jj) * R(abs(jj - ii), m) + H3_2(1, ii, jj) / (2 * m);
        A(3 + ii + 2 * m, jj + 2 * m) = H3(1, 2, ii, jj) / (2 * m);
        A(3 + ii + 2 * m, jj + 4 * m) = H4_1(1, ii, jj) * R(abs(jj - ii), m) + H4_2(1, ii, jj) / (2 * m);
        A(3 + ii + 2 * m, jj + 6 * m) = H4(1, 2, ii, jj) / (2 * m);
   
        % 3d equation kernels
        A(3 + ii + 4 * m, jj) = H3(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 4 * m, jj + 2 * m) = H3_1(2, ii, jj) * R(abs(jj - ii), m) + H3_2(2, ii, jj) / (2 * m);
        A(3 + ii + 4 * m, jj + 4 * m) = H4(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 4 * m, jj + 6 * m) = H4_1(2, ii, jj) * R(abs(jj - ii), m) + H4_2(2, ii, jj) / (2 * m);
        
        % 4th equation kernels
        A(3 + ii + 6 * m, jj) = H5(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 6 * m, jj + 2 * m) = H5_1(2, ii, jj) * R(abs(jj - ii), m) + H5_2(2, ii, jj) / (2 * m);
        A(3 + ii + 6 * m, jj + 4 * m) = H6(2, 1, ii, jj) / (2 * m);
        A(3 + ii + 6 * m, jj + 6 * m) = H6(2, 2, ii, jj) / (2 * m);
    end
    
    %---------------------------%
end
    
%---------------------------%

% Write data to file 

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

%---------------------------%

%disp('the end of matrix initialization');
%main