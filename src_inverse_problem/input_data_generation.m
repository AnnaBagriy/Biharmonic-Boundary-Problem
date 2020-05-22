format long

tic

% Constants initialization

m = 16;

A0 = 1;
A1 = 1;
A2 = 1;

nu = 0.5;

y1_point = 0;
y2_point = 0;

x1_point = 1.5;
x2_point = 1;

r_2 = @(x1, y1, x2, y2) sqrt((x1 - y1).^2 + (x2 - y2).^2);

% g2 = @(x, y) x + y;
% q2 = @(x, y) x - y;

disp([newline char(9) 'RESULTS'])
disp(['m = ' num2str(m)])

% Initialize system of linear equations

dps = DirectProblemSolver(m, A0, A1, A2, nu);

% Set exact domain

%r_ex = @(t) 1 + 0.15 * cos(3 * t);
r_ex = @(t) sqrt(cos(t).^2 + 0.25 * sin(t).^2);
%r = @(t) 0.9;

x1_ex = @(t) r_ex(t) .* cos(t);
y1_ex = @(t) r_ex(t) .* sin(t);

r_2_ex = @(t) ((cos(t) ./ 2).^10 + (2 * sin(t) / 3).^10).^-0.1;

x2 = @(t) r_2_ex(t) .* cos(t);
y2 = @(t) r_2_ex(t) .* sin(t);

dps.SetBoundary1(x1_ex, y1_ex);
dps.SetBoundary2(x2, y2);

% Right vector functions
% g2 = @ (x1, x2) ((dps.n1_on_2 .* (x1 - y1_point) + dps.n1_on_2 .* (x2 - y2_point)) .* (1 + 2 * log(r_2(x1, y1_point, x2, y2_point)))) / (8 * pi);
% 
% q2 = @ (x1, x2) nu * 6 + 2 * log(r_2(x1, y1_point, x2, y2_point).^2) + (1 - nu) *...
%     (dps.n1_on_2.^2 .* (1 + 2 * (x1 - y1_point).^2 ./ (r_2(x1, y1_point, x2, y2_point).^2) + log(r_2(x1, y1_point, x2, y2_point).^2)) +...
%      dps.n2_on_2.^2 .* (1 + 2 * (x2 - y2_point).^2 ./ (r_2(x1, y1_point, x2, y2_point).^2) + log(r_2(x1, y1_point, x2, y2_point).^2)) +...
%      4 * dps.n1_on_2 .* dps.n2_on_2 .* (x1 - y1_point) .* (x2 - y2_point) ./ (r_2(x1, y1_point, x2, y2_point).^2)) * 0.1;

g2 = @(x1, x2) 2 * x1 .* dps.n1_on_2 - 2 * x2 .* dps.n2_on_2;
q2 = @(x1, x2) 2 * (1 - nu) * (dps.n1_on_2.^2 - dps.n2_on_2.^2);
% 
%  g2 = @(x1, x2) dps.n1_on_2 + dps.n2_on_2;
%  q2 = @(x1, x2) x1 * 0;
% 
%  g2 = @(x, y) x - y;
%  q2 = @(x, y) x + y;

dps.SetBoundaryFunctions(g2, q2);

% Solve system of the equations
% Find densities and f on Г2

dps.FindDensitiesAndConstants();
f = dps.FindF();

disp([newline 'f(x) (x on Г2) = ' num2str(f) newline]);

u_ex = dps.FindExactU(x1_point, x2_point, y1_point, y2_point);
%u_ex = -1.9519;
u = dps.FindU(x1_point, x2_point);

error = u_ex - u;

disp([newline 'u_ex(x) = ' num2str(u_ex)]);
disp([newline 'u_ap(x) = ' num2str(u)]);
disp([newline 'error   = ' num2str(error)]);
% Plot domain

tn = 0:0.01:2 * pi;
% 
%  h = figure;
% % 
plot(x1_ex(tn), y1_ex(tn), 'k');
hold on;
plot(x2(tn), y2(tn), 'k');
hold on;
% 
% saveas(h,'sample3.pdf')
 
toc
