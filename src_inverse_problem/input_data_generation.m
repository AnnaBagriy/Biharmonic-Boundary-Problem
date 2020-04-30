
tic

% Constants initialization

m = 32;

A0 = 1;
A1 = 1;
A2 = 1;

nu = 0.5;

% Right vector functions
g2 = @(x, y) x + y;
q2 = @(x, y) x - y;

disp([newline char(9) 'RESULTS'])
disp(['m = ' num2str(m)])

% Initialize system of linear equations

mx = DirectProblemSolver(m, g2, q2, A0, A1, A2, nu);

% Set exact domain

r = @(t) sqrt(cos(t).^2 + 0.25 * sin(t).^2);

x1 = @(t) r(t) .* cos(t);
y1 = @(t) r(t) .* sin(t);

x2 = @(t) 2 * cos(t);
y2 = @(t) 2 * sin(t);

mx.SetBoundary1(x1, y1, 0);
mx.SetBoundary2(x2, y2);

% Solve system of the equations
% Find densities and f on Г2

mx.FindDensitiesAndConstants();
mx.FindF();

% Plot domain

tn = 0:0.01:2 * pi;

plot(x1(tn), y1(tn));
hold on;
plot(x2(tn), y2(tn));
hold on;

toc
