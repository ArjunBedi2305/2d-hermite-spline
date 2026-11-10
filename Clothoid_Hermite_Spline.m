clear; clc; close all;

x0 = 0;   y0 = 0;   theta0 = deg2rad(90);
x1 = 12;  y1 = 3;   theta1 = deg2rad(135);

d     = hypot(x1 - x0, y1 - y0);
phi   = atan2(y1 - y0, x1 - x0);
wrap  = @(a) mod(a + pi, 2*pi) - pi;
phi0  = wrap(theta0 - phi);
phi1  = wrap(theta1 - phi);
delta = wrap(phi1 - phi0);

th = @(u, A) A*u.^2 + (delta - A)*u + phi0;
g  = @(A) integral(@(u) sin(th(u, A)), 0, 1);
hf = @(A) integral(@(u) cos(th(u, A)), 0, 1);

% Scan for sign changes of g, keep roots with positive length, pick the shortest
As = linspace(-60, 60, 1201);
gs = arrayfun(g, As);
best = NaN; bestLen = inf;
for i = find(sign(gs(1:end-1)) ~= sign(gs(2:end)))
    Ai = fzero(g, [As(i) As(i+1)]);
    hi = hf(Ai);
    if hi > 1e-6 && d/hi < bestLen
        best = Ai; bestLen = d/hi;
    end
end
if isnan(best), error('No clothoid solution found for these poses.'); end

A    = best;
Ltot = bestLen;
k0   = (delta - A) / Ltot;
dk   = 2*A / Ltot^2;

N = 1000;
s = linspace(0, Ltot, N);
heading = theta0 + k0*s + 0.5*dk*s.^2;
x_s = x0 + cumtrapz(s, cos(heading));
y_s = y0 + cumtrapz(s, sin(heading));
kappa = k0 + dk*s;

fprintf('Arc length            = %.3f\n', Ltot);
fprintf('Initial curvature     = %.4f\n', kappa(1));
fprintf('Final curvature       = %.4f\n', kappa(end));
fprintf('Curvature rate dk/ds  = %.4f\n', dk);
fprintf('Peak |curvature|      = %.4f\n', max(abs(kappa)));
fprintf('Min turn radius       = %.3f\n', 1/max(abs(kappa)));
fprintf('End position error    = %.2e\n', hypot(x_s(end) - x1, y_s(end) - y1));
fprintf('End heading error     = %.2e rad\n', wrap(heading(end) - theta1));

figure('Color', 'w');
hold on; grid on; box on;
plot(x_s, y_s, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Clothoid Path');
plot(x0, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Start');
plot(x1, y1, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'End');
a = 0.15 * d;
quiver(x0, y0, a*cos(theta0), a*sin(theta0), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
quiver(x1, y1, a*cos(theta1), a*sin(theta1), 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
xlabel('X'); ylabel('Y');
title('Clothoid Path');
legend('Location', 'best');
axis equal;
