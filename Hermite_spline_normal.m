%% 2D Cubic Hermite Spline Parameterization with Shared Tangent Length L
clear; clc; close all;

%Define Input Poses (Positions and Angles in Radians)
x0 = 0;
y0 = 0;
theta0 = deg2rad(90);   % Straight up (North) % Define heading of initial vector in degrees here 

x1 = 12;
y1 = 3;
theta1 = deg2rad(135); % Define heading of final vector in degrees here 

% Single shared tangent magnitude (velocity factor)
L = 5.0; 

% Pre-calculate Trigonometric Components
cos0 = cos(theta0);
sin0 = sin(theta0);
cos1 = cos(theta1);
sin1 = sin(theta1);

%Compute Polynomial Coefficients (ax*t^3 + bx*t^2 + cx*t + dx)
% X-axis coefficients
ax = 2*(x0 - x1) + L*(cos0 + cos1);
bx = 3*(x1 - x0) - L*(2*cos0 + cos1);
cx = L*cos0;
dx = x0;

% Y-axis coefficients
ay = 2*(y0 - y1) + L*(sin0 + sin1);
by = 3*(y1 - y0) - L*(2*sin0 + sin1);
cy = L*sin0;
dy = y0;

%Evaluate the Spline Curve over t = [0, 1]
t = linspace(0, 1, 200); % 200 points for a smooth curve

x_t = ax.*t.^3 + bx.*t.^2 + cx.*t + dx;
y_t = ay.*t.^3 + by.*t.^2 + cy.*t + dy;

%Plot the Results
figure('Color', 'w');
hold on; grid on; box on;

% Plot the generated spline curve
plot(x_t, y_t, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Hermite Spline Path');

% Plot start and end control points
plot(x0, y0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'P_0 (Start)');
plot(x1, y1, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'P_1 (End)');

% Draw tangent vectors as arrows to visualize orientation
quiver(x0, y0, L*cos0, L*sin0, 0, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.1, 'DisplayName', 'm_0 (Start Tangent)');
quiver(x1, y1, L*cos1, L*sin1, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.1, 'DisplayName', 'm_1 (End Tangent)');

%Formatting the Plot
xlabel('X Coordinate', 'FontSize', 12);
ylabel('Y Coordinate', 'FontSize', 12);
title('2D Cubic Hermite Spline Simulation', 'FontSize', 14);
legend('Location', 'best');
axis equal; % Ensures grid units are square and unskewed

% Extend axes slightly past data bounds for visibility
xlim([min(x_t)-3, max(x_t)+5]);
ylim([min(y_t)-3, max(y_t)+5]);
