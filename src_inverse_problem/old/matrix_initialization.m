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

%---------------------------%

% m = 4;
% h = pi / m;
% jj = 1:2 * m;
% s = (jj - 1) * h;
% 
% nu = 0.5;

x1 = x_vector(1, s, size(s, 2));
y1 = y_vector(1, s, size(s, 2));
x2 = x_vector(2, s, size(s, 2));
y2 = y_vector(2, s, size(s, 2));

n_11 = n1(1, s);
n_21 = n2(1, s);
n_12 = n1(2, s);
n_22 = n2(2, s);

A0 = 1;
A1 = 1;
A2 = 1;

% Right vector functions
g2 = @(x, y) x + y;
q2 = @(x, y) x - y;

%---------------------------%

% Derivatives

syms t;

y_derivative1_sym = matlabFunction(diff(y_vector(1, t, size(t)), t));
y_derivative2_sym = matlabFunction(diff(y_vector(2, t, size(t)), t));

y_2derivative1_sym = matlabFunction(diff(diff(y_vector(1, t, size(t)), t), t));
y_2derivative2_sym = matlabFunction(diff(diff(y_vector(2, t, size(t)), t), t));

x_derivative1_sym = matlabFunction(diff(x_vector(1, t, size(t)), t));
x_derivative2_sym = matlabFunction(diff(x_vector(2, t, size(t)), t));

x_2derivative1_sym = matlabFunction(diff(diff(x_vector(1, t, size(t)), t), t));
x_2derivative2_sym = matlabFunction(diff(diff(x_vector(2, t, size(t)), t), t));

y_derivative1 = y_derivative1_sym(s);
y_derivative2 = y_derivative2_sym(s);

x_derivative1 = x_derivative1_sym(s);
x_derivative2 = x_derivative2_sym(s);

x_2derivative1 = x_2derivative1_sym(s);
x_2derivative2 = x_2derivative2_sym(s);

y_2derivative1 = y_2derivative1_sym(s);
y_2derivative2 = y_2derivative2_sym(s);


% mx = MainSystem(m, g2, q2, A0, A1, A2, nu);
% 
% A = mx.InitializeMatrix();
% y = mx.InitializeRightVector();
