
% Constants initialization

m = 8;

A0 = 1;
A1 = 1;
A2 = 1;

nu = 0.5;

epsilon = 1e-3;

% Right vector functions
g2 = @(x, y) x + y;
q2 = @(x, y) x - y;

disp([newline char(9) 'RESULTS'])
disp(['m = ' num2str(m)])

% Initialize system of linear equations

dps = DirectProblemSolver(m, g2, q2, A0, A1, A2, nu);

% Set initial domain for approximation

 % Approximate Г1
r = @(t) 1;

x1 = @(t) r(t) .* cos(t);
y1 = @(t) r(t) .* sin(t);

 % Exact Г2
x2 = @(t) 2 * cos(t);
y2 = @(t) 2 * sin(t);

dps.SetBoundary1(x1, y1, 0);
dps.SetBoundary2(x2, y2);

% Plot approximate domain

tn = 0:0.01:2 * pi;

plot(x1(tn), y1(tn));
hold on;

% Solve system of the equations
% Find densities and f on Г2

[fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2] = dps.FindDensitiesAndConstants();
f = dps.FindF();

% Solve incorrect equation

 % Constants initialization
n = 4;
lambda = 1e-10;

 % Initialization
ies = IncorrectEquationSolver(n, m, lambda, f);

ies.SetBoundaryAndNormal(dps.s, dps.x1, dps.y1, dps.x2, dps.y2, dps.n1_on_1, dps.n2_on_1, dps.n1_on_2, dps.n2_on_2, 0);
ies.InitializeDensitiesAndConstants(fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2);

 % Solving
q_m = ies.Getq_m();
    
 % Results displaying
disp(['||q_m|| = ' num2str(norm(q_m))]);
    
while abs(norm(q_m)) >= epsilon
    
    % Set initial domain for approximation

     % Approximate Г1
    r = @(t) r(t);

    x1 = @(t) r(t) .* cos(t);
    y1 = @(t) r(t) .* sin(t);

    dps.SetBoundary1(x1, y1, q_m);
    
    % Plot approximate domain

%     plot(x1(dps.s), y1(dps.s));
%     hold on;

    % Solve system of the equations
    % Find densities and f on Г2

    [fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2] = dps.FindDensitiesAndConstants();
    f = dps.FindF();

    % Solve incorrect equation

     % Constants initialization
    n = 4;
    lambda = 1e-10;

     % Initialization
    ies = IncorrectEquationSolver(n, m, lambda, f);

    ies.SetBoundaryAndNormal(dps.s, dps.x1, dps.y1, dps.x2, dps.y2, dps.n1_on_1, dps.n2_on_1, dps.n1_on_2, dps.n2_on_2, q_m);
    ies.InitializeDensitiesAndConstants(fi_func_1, fi_func_2, psi_func_1, psi_func_2, a0, a1, a2);

     % Solving
    q_m = ies.Getq_m();

     % Results displaying
    disp(['||q_m|| = ' num2str(norm(q_m))]);
    
end
