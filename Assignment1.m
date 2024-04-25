%% Advanced Astrodynamics (ASEN6060) - Spring 2024
% Instructor:   Prof. Bosanac
% Student:      Giovanni Fereoli
% Student ID:   111040314

%% Problem 1
clear; clc; close;

% Data
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
GM_sun = 132712440041.93938;
G_tilde = 6.67408*1e-20; % km^3/kg*s^2
e_earth = 0.016708617; 
e_moon = 0.05490; 
a_earth = 149598023; % km, wrt sun
a_moon = 384400;  % km, wrt moon

% Computation m^* in kg
M_tilde_earth = GM_earth / G_tilde;
M_tilde_moon = GM_moon / G_tilde;
M_tilde_sun = GM_sun / G_tilde;
m_star_se = M_tilde_sun + M_tilde_earth; 
m_star_em = M_tilde_earth + M_tilde_moon;

% Computation l^* in km
l_star_se = a_earth; % OSS: approximation
l_star_em = a_moon;

% Computation t^* in seconds
t_star_se = sqrt(l_star_se^3 / (G_tilde * m_star_se));
t_star_em = sqrt(l_star_em^3 / (G_tilde * m_star_em));

% Computations mass ratios
mu_se = GM_earth / (GM_sun + GM_earth); % Mass Ratio
mu_em = GM_moon / (GM_earth + GM_moon);

%% Problem 2
clc; close;

%Initial conditions
x0_1 = [0.98; 0; 0; 0; 1.2; 0];
x0_2 = [0.98; 0; 0; 0; 1.7; 0];
x0_3 = [0.12; 0; 0; 0; 3.45; 0];
x0_4 = [0.12; 0; 0; 0; 3.48; 0];

% Time of flights
ToF_1 = 2;
ToF_2 = 8;
ToF_3 = 25;
ToF_4 = 25;

% Initialization
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
mu = GM_moon / (GM_earth + GM_moon);
C1 = JacobiConstant(x0_1, mu);
C2 = JacobiConstant(x0_2, mu);
C3 = JacobiConstant(x0_3, mu);
C4 = JacobiConstant(x0_4, mu);
options = odeset('RelTol', 2.22045*1e-14, 'AbsTol', 2.22045*1e-16);

% Integration
[~, xx1] = ode113(@(t,X) CRTBP(t, X, mu), [0 ToF_1], x0_1, options);
[~, xx2] = ode113(@(t,X) CRTBP(t, X, mu), [0 ToF_2], x0_2, options);
[~, xx3] = ode113(@(t,X) CRTBP(t, X, mu), [0 ToF_3], x0_3, options);
[~, xx4] = ode113(@(t,X) CRTBP(t, X, mu), [0 ToF_4], x0_4, options);

% Plots
gca = figure(1);
subplot(2,2,1);
plot3(xx1(:,1), xx1(:,2), xx1(:,3), 'b', 'LineWidth', 1);
hold on;
ZVCxy(mu, C1);
plot3(1-mu, 0, 0, 'c.', 'MarkerSize', 10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
title(['First $\mathbf{x}_{0}$ with $C_J$=', num2str(C1)],...
    'interpreter','latex');
legend('Trajectory', 'ZVC', 'Moon', 'interpreter','latex',...
    'Location','northeast', 'FontSize', 5);
grid on;
real axis;
xlim([0.97, 1.005]);
ylim([- 0.05, 0.05]);
view(2);
subplot(2,2,2);
plot3(xx2(:,1), xx2(:,2), xx2(:,3), 'r', 'LineWidth', 1);
hold on;
ZVCxy(mu, C2);
plot3(1-mu, 0, 0, 'c.', 'MarkerSize', 10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
title(['Second $\mathbf{x}_{0}$ with $C_J$=', num2str(C2)],...
    'interpreter','latex');
legend('Trajectory','ZVC', 'Moon', 'interpreter','latex',...
    'Location','northeast', 'FontSize', 5);
hold on;
grid on;
real axis;
xlim([0.8, 1.2]);
ylim([- 0.2, 0.2]);
view(2);
subplot(2,2,3);
plot3(xx3(:,1), xx3(:,2), xx3(:,3), 'g', 'LineWidth', 1);
hold on;
ZVCxy(mu, C3);
plot3(-mu, 0, 0, 'y.', 'MarkerSize', 10);
plot3(1-mu, 0, 0, 'c.', 'MarkerSize', 10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
title(['Third $\mathbf{x}_{0}$ with $C_J$=', num2str(C3)],...
    'interpreter','latex');
legend('Trajectory','ZVC', 'Earth', 'Moon', 'interpreter','latex',...
    'Location','northeast', 'FontSize', 5);
grid on;
real axis;
view(2);
subplot(2,2,4);
plot3(xx4(:,1), xx4(:,2), xx4(:,3), 'm', 'LineWidth', 1);
hold on;
ZVCxy(mu, C4);
plot3(-mu, 0, 0, 'y.', 'MarkerSize', 10);
plot3(1-mu, 0, 0, 'c.', 'MarkerSize', 10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
title(['Fourth $\mathbf{x}_{0}$ with $C_J$=', num2str(C4)],...
    'interpreter','latex');
legend('Trajectory', 'ZVC', 'Earth', 'Moon', 'interpreter','latex',...
    'Location','northeast', 'FontSize', 5);
grid on;
real axis;
view(2);

% Save plot
exportgraphics(gca, 'Prob2.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Last computations
x0_tilde_3 = zeros(6);
x0_tilde_3(1:3) = x0_3(1:3) * l_star_em;
x0_tilde_3(4:6) = x0_3(4:6) * l_star_em / t_star_em;
ToF_tilde_3 = ToF_3 * t_star_em;
earth_pos = [-mu; 0; 0];
rho_earth_x0_3 = norm(earth_pos - x0_3(1:3)) * l_star_em;
periods_ToF_3 = ToF_3 / (2*pi);

%% Problem 2C
clear; clc; close;

%Initial conditions
x0_1 = [0.98; 0; 0; 0; 1.2; 0];
x0_2 = [0.98; 0; 0; 0; 1.7; 0];
x0_3 = [0.12; 0; 0; 0; 3.45; 0];
x0_4 = [0.12; 0; 0; 0; 3.48; 0];

% Time of flights
ToF_1 = 2;
ToF_2 = 8;
ToF_3 = 25;
ToF_4 = 25;

% Initialization
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
mu = GM_moon / (GM_earth + GM_moon);
options = odeset('RelTol', 2.22045*1e-14, 'AbsTol', 2.22045*1e-16);

% Integration
len = 10000;
[t1, xx1] = ode113(@(t,X) CRTBP(t, X, mu), ...
    linspace(0, ToF_1, len), x0_1, options);
[t2, xx2] = ode113(@(t,X) CRTBP(t, X, mu), ...
    linspace(0, ToF_2, len), x0_2, options);
[t3, xx3] = ode113(@(t,X) CRTBP(t, X, mu),...
    linspace(0, ToF_3, len), x0_3, options);
[t4, xx4] = ode113(@(t,X) CRTBP(t, X, mu), ...
    linspace(0, ToF_4, len), x0_4, options);

% Jacobi Constants
C1_vec = zeros(len, 1);
C2_vec = zeros(len, 1);
C3_vec = zeros(len, 1);
C4_vec = zeros(len, 1);
cumRelStd_C1 = zeros(len, 1);
cumRelStd_C2 = zeros(len, 1);
cumRelStd_C3 = zeros(len, 1);
cumRelStd_C4 = zeros(len, 1);
for i = 1:len
    % First initial condition
    C1_vec(i) = JacobiConstant(xx1(i,:), mu);
    cumRelStd_C1(i) = std(C1_vec(1:i)) / mean(C1_vec(1:i));

    % Second initial condition
    C2_vec(i) = JacobiConstant(xx2(i,:), mu);
    cumRelStd_C2(i) = std(C2_vec(1:i)) / mean(C2_vec(1:i));

    % Third initial condition
    C3_vec(i) = JacobiConstant(xx3(i,:), mu);
    cumRelStd_C3(i) = std(C3_vec(1:i)) / mean(C3_vec(1:i));

    % Fourth initial condition
    C4_vec(i) = JacobiConstant(xx4(i,:), mu);
    cumRelStd_C4(i) = std(C4_vec(1:i)) / mean(C4_vec(1:i));
end

%Plot
samples = 1:len;
figure(1)
semilogy(samples, cumRelStd_C1, 'm', 'LineWidth', 1);
hold on;
semilogy(samples, cumRelStd_C2, 'r', 'LineWidth', 1);
semilogy(samples, cumRelStd_C3, 'b', 'LineWidth', 1);
semilogy(samples, cumRelStd_C4, 'g', 'LineWidth', 1);
xlabel('Samples [-]', 'interpreter','latex');
ylabel('$\sigma_{C}/\mu_{C}$ [-]', 'interpreter','latex');
legend('First $\mathbf{x}_0$, Integration time: 2 [-]',...
    'Second $\mathbf{x}_0$, Integration time: 8 [-]',...
    'Third $\mathbf{x}_0$, Integration time: 25 [-]',...
    'Fourth $\mathbf{x}_0$, Integration time: 25 [-]',...
    'interpreter','latex', 'Location','southeast', 'FontSize', 9);
grid on;
real axis;

% Save plot
exportgraphics(gca, 'Prob2C.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 3
clear; clc; close;

% Initial condition and Time of Flight
x0_3 = [0.12; 0; 0; 0; 3.45; 0];
ToF_3 = 25;

% Initialization
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
mu = GM_moon / (GM_earth + GM_moon);
C3 = JacobiConstant(x0_3, mu);
options = odeset('RelTol', 2.22045*1e-14, 'AbsTol', 2.22045*1e-16,...
        'Events',@event);

% Integration 
[tt3, xx3] = ode113(@(t,X) CRTBP(t, X, mu), [0 ToF_3], x0_3, options);

% Plots
gca = figure(1);
plot3(xx3(:,1), xx3(:,2), xx3(:,3), 'g', 'LineWidth', 1);
hold on;
ZVCxy(mu, C3);
plot3(-mu, 0, 0, 'y.', 'MarkerSize', 10);
plot3(1-mu, 0, 0, 'c.', 'MarkerSize', 10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
title('Third $\mathbf{x}_0$ with event $y>0, \ \dot{y}>0$',...
    'interpreter','latex');
legend('Trajectory', 'ZVC', 'Earth', 'Moon', 'interpreter','latex',...
    'Location','northeast');
view(2)
grid on;
real axis;

% Save plot
exportgraphics(gca, 'Prob3.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Results at event
xx3_event = xx3(end, :);
tt3_event = tt3(end);

%% Problem 4
clear; clc; close;

% Initialization
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
mu = GM_moon / (GM_earth + GM_moon);

% Jacobi Constant
C1 = 3.189;
C2 = 3.173;
C3 = 3.013;
C4 = 2.995;

% Plots
gca = figure(1);
subplot(2,2,1);
ZVCxy(mu, C1);
hold on;
plot(-mu, 0, 'b.','MarkerSize',10);
plot(1-mu, 0, 'r.','MarkerSize',10);
legend('ZVC', 'Earth', 'Moon', 'interpreter','latex',...
    'Location','northeast', 'FontSize', 7);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$C_J = 3.189$', 'interpreter','latex');
grid on;
real axis;
subplot(2,2,2);
ZVCxy(mu, C2);
hold on;
plot(-mu, 0, 'b.','MarkerSize',10);
plot(1-mu, 0, 'r.','MarkerSize',10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$C_J = 3.173$', 'interpreter','latex');
grid on;
real axis;
subplot(2,2,3);
ZVCxy(mu, C3);
hold on;
plot(-mu, 0, 'b.','MarkerSize',10);
plot(1-mu, 0, 'r.','MarkerSize',10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$C_J = 3.013$', 'interpreter','latex');
grid on;
real axis;
subplot(2,2,4);
ZVCxy(mu, C4);
hold on;
plot(-mu, 0, 'b.','MarkerSize',10);
plot(1-mu, 0, 'r.','MarkerSize',10);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$C_J = 2.995$', 'interpreter','latex');
grid on;
real axis;

% Save plot
exportgraphics(gca, 'Prob4.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 4D
clear; clc; close;

% Initialization
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
mu = GM_moon / (GM_earth + GM_moon);

% Jacobi Constant
C1 = 3.18841;
C2 = 3.1722;
C3 = 3.0124;
C45 = 2.988;

% Plots
gca = figure(1);
subplot(2,2,1);
ZVCxy(mu, C1);
hold on;
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$L_1$ Position', 'interpreter','latex');
grid on;
real axis;
xlim([0.832,0.842]);
ylim([-0.02,0.02]);
subplot(2,2,2);
ZVCxy(mu, C2);
hold on;
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$L_2$ Position', 'interpreter','latex');
grid on;
real axis;
xlim([1.151,1.161]);
ylim([-0.02,0.02]);
subplot(2,2,3);
ZVCxy(mu, C3);
hold on;
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$L_3$ Position', 'interpreter','latex');
grid on;
real axis;
xlim([-1.05,-0.95]);
ylim([-0.01,0.01]);
subplot(2,2,4);
ZVCxy(mu, C45);
hold on;
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
title('$L_4$ Position', 'interpreter','latex');
grid on;
real axis;
xlim([0.44,0.54]);
ylim([0.861,0.871]);

% Know position Lagrangian points from ZVC plot
[x_L, y_L] = ginput(4);

% Add to plot Lgrange Points
figure(1);
subplot(2,2,1);
plot(x_L(1), y_L(1), 'r.','MarkerSize',20);
legend('ZVC', 'Lagrange Point', 'interpreter','latex',...
    'Location','northeast', 'FontSize', 7);
subplot(2,2,2);
plot(x_L(2), y_L(2), 'r.','MarkerSize',20);
subplot(2,2,3);
plot(x_L(3), y_L(3), 'r.','MarkerSize',20);
subplot(2,2,4);
plot(x_L(4), y_L(4), 'r.','MarkerSize',20);

% Save plot
exportgraphics(gca, 'Prob4D.pdf', 'ContentType','image',...
     'Resolution', 1000);

%% Functions 

% Compute ZVC in the xy-plane
function ZVCxy(mu, C)

% Initialization 
x_zvc = -1.5:0.001:1.55;
y_zvc = -1.5:0.001:1.55;
Z_zvc = zeros(length(y_zvc), length(x_zvc));

% ZVS computations
for i=1:length(x_zvc)
    for j=1:length(y_zvc)
        Z_zvc(j, i)= (x_zvc(i)^2+y_zvc(j)^2) + ...
            2*(1-mu)./sqrt((x_zvc(i)+mu)^2+y_zvc(j)^2) + ...
            2*mu./sqrt((x_zvc(i)-1+mu)^2+y_zvc(j)^2);
    end
end

% Plot
contourf(x_zvc, y_zvc, -Z_zvc,[-C -C]);

end

% Compute Jacobi Constant
function C = JacobiConstant(X, mu)

% Initialization
x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);

% Jacobi Constant Computation
C = (x^2+y^2) + 2*(1-mu)/sqrt((x+mu)^2+y^2+z^2) + ...
    2*mu/sqrt((x-1+mu)^2+y^2+z^2) - sqrt(xdot^2+ydot^2+zdot^2)^2;

end

% CR3BP Equations of Motions
function dXdt = CRTBP(~, X, mu)

%Initialize
dXdt = zeros(6,1);

x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);

% CRTBP dynamics
r1_norm = sqrt((x+mu)^2+y^2+z^2);
r2_norm = sqrt((x+mu-1)^2+y^2+z^2);

dXdt(1:3) = [xdot; ydot; zdot];
dXdt(4:6) = [2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3;...
    -(1-mu)*z/r1_norm^3-mu*z/r2_norm^3];
end

% Event
function [value, isterminal, direction] = event(~, X, ~)

value = X(2);   % y = 0
isterminal = 1;
direction = 1; % ydot > 0, event function increasing

end