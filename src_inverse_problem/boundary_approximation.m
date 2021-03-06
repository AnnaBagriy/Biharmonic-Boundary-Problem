
tic

% Generate f on Г2

input_data_generation

% Constants initialization

epsilon = 1e-3;

n = 4;
lambda = 1e-3;

disp([newline char(9) 'RESULTS'])
disp(['m = ' num2str(m) newline])

% Initialize direct problem solver

dps = DirectProblemSolver(m, g2, q2, A0, A1, A2, nu);

% Set initial domain for approximation

 % Approximate Г1
r = @(t) 1.5;

x1 = @(t) cos(t);
y1 = @(t) sin(t);

x1_r = @(t) r(t) .* cos(t);
y1_r = @(t) r(t) .* sin(t);

dps.SetBoundary1(x1_r, y1_r);
dps.SetBoundary2(x2, y2);

% Plot approximate domain

% tn = 0:0.01:2 * pi;
% 
plot(x1_r(tn), y1_r(tn), 'r--');
hold on;

% Solve system of the equations
% Find densities and f on Г2

[fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2] = dps.FindDensitiesAndConstants();

% Solve incorrect equation

 % Initialization
ies = IncorrectEquationSolver(n, m, lambda, f);

ies.SetBoundaryAndNormal(dps.s, dps.x1, dps.y1, dps.x2, dps.y2, dps.n1_on_1, dps.n2_on_1, dps.n1_on_2, dps.n2_on_2);
ies.SetDerivatives(dps.x_derivative_1, dps.y_derivative_1, dps.x_derivative_2, dps.y_derivative_2, dps.x_second_derivative_1,...
                   dps.y_second_derivative_1, dps.x_second_derivative_2, dps.y_second_derivative_2);

ies.InitializeDensitiesAndConstants(fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2);

 % Solving
q_m = ies.Getq_m();
    
 % Results displaying
disp(['q_m = ' num2str(q_m)]);
disp(['||q_m|| = ' num2str(norm(q_m))]);

r_q_m = r(dps.s);

ii = 1;
% while abs(norm(q_m)) >= epsilon
% 
%     % Approximate boundary by applying r = r + q
%     r_q_m = dps.SetApproximationBoundary1(x1, y1, r_q_m, q_m);
% 
%     % Solve system of the equations
%     % Find densities and f on Г2
% 
%     [fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2] = dps.FindDensitiesAndConstants();
% 
%     % Solve incorrect equation
% 
%     ies.SetBoundaryAndNormal(dps.s, dps.x1, dps.y1, dps.x2, dps.y2, dps.n1_on_1, dps.n2_on_1, dps.n1_on_2, dps.n2_on_2);
%     ies.SetDerivatives(dps.x_derivative_1, dps.y_derivative_1, dps.x_derivative_2, dps.y_derivative_2, dps.x_second_derivative_1, dps.y_second_derivative_1, dps.x_second_derivative_2, dps.y_second_derivative_2);
% 
%     ies.InitializeDensitiesAndConstants(fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2);
% 
%      % Solving
%     q_m = ies.Getq_m();
% 
%      % Results displaying
% %      disp([newline 'x1 = ' num2str(dps.x1)]);
% %      disp(['y1 = ' num2str(dps.y1)]);
%        disp([newline 'q_m = ' num2str(q_m)]);
%        disp(['||q_m|| = ' num2str(norm(q_m))]);
% 
%     ii = ii + 1;
% end

disp([newline 'x_ex     = ', num2str(x1_ex(dps.s))])
disp(['x_approx = ', num2str(dps.x1)])

disp([newline 'y_ex     = ', num2str(y1_ex(dps.s))])
disp(['y_approx = ', num2str(dps.y1)])

x_error = norm(x1_ex(dps.s) - dps.x1);
y_error = norm(y1_ex(dps.s) - dps.y1);

disp([newline 'x_error = ', num2str(x_error)])
disp(['y_error = ', num2str(y_error)])

plot(dps.x1, dps.y1, 'b--');
hold on;

disp([newline 'Iterations = ', num2str(ii)]);

toc
