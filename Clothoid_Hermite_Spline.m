%% 2D Clothoid G1 Interpolation Between Two Poses
clear; clc; close all;

%Define Input Poses (Positions and Angles in Radians)
x0 = 0;   y0 = 0;   theta0 = deg2rad(90);
x1 = 12;  y1 = 3;   theta1 = deg2rad(135);

%Normalize the Problem
dx = x1 - x0;
dy = y1 - y0;
d   = hypot(dx, dy);          % chord length
phi = atan2(dy, dx);          % chord direction

wrap = @(a) mod(a + pi, 2*pi) - pi;   % wrap angle to [-pi, pi)
phi0 = wrap(theta0 - phi);    % start heading relative to chord
phi1 = wrap(theta1 - phi);    % end heading relative to chord
delta = wrap(phi1 - phi0);    % total turning angle

%Solve for the Clothoid Parameter A
% Normalized heading over u in [0,1]: theta(u) = A*u^2 + (delta - A)*u + phi0
% The curve reaches the target point when the integral of sin(theta) is 0.
th = @(u, A) A*u.^2 + (delta - A)*u + phi0;
g  = @(A) integral(@(u) sin(th(u, A)), 0, 1);

A0 = 3*(phi0 + phi1);         % good initial guess (exact when phi0+phi1 = 0)
A  = fzero(g, A0);

% Recover Length and Curvature Parameters
h     = integral(@(u) cos(th(u, A)), 0, 1);
Ltot  = d / h;                % total arc length
k0    = (delta - A) / Ltot;   % initial curvature
dk    = 2*A / Ltot^2;         % curvature rate (clothoid sharpness)

% Evaluate the Curve Over Arc Length s
N = 500;
s = linspace(0, Ltot, N);

heading = theta0 + k0*s + 0.5*dk*s.^2;
x_s = x0 + cumtrapz(s, cos(heading));
y_s = y0 + cumtrapz(s, sin(heading));
kappa = k0 + dk*s;            % curvature varies linearly with s

fprintf('Arc length = %.3f, k0 = %.4f, k1 = %.4f\n', Ltot, kappa(1), kappa(end));
fprintf('End point error = %.2e\n', hypot(x_s(end)-x1, y_s(end)-y1));
fprintf('End heading error = %.2e rad\n', wrap(heading(end) - theta1));

%Plot the Path
figure('Color', 'w');
hold on; grid on; box on;
plot(x_s, y_s, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Clothoid Path');
plot(x0, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'P_0 (Start)');
plot(x1, y1, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'P_1 (End)');

arrowLen = 0.15 * d;
quiver(x0, y0, arrowLen*cos(theta0), arrowLen*sin(theta0), 0, 'k', ...
       'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'Start Heading');
quiver(x1, y1, arrowLen*cos(theta1), arrowLen*sin(theta1), 0, 'g', ...
       'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'End Heading');

xlabel('X Coordinate', 'FontSize', 12);
ylabel('Y Coordinate', 'FontSize', 12);
title('2D Clothoid G1 Interpolation', 'FontSize', 14);
legend('Location', 'best');
axis equal;
xlim([min(x_s)-1, max(x_s)+2]);
ylim([min(y_s)-1, max(y_s)+2]);

% Plot Curvature vs Arc Length
figure('Color', 'w');
plot(s, kappa, 'r-', 'LineWidth', 2); grid on;
xlabel('Arc length s'); ylabel('Curvature \kappa (1/unit)');
title('Curvature Along Path (linear in s)');
