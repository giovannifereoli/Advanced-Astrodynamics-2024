%% Advanced Astrodynamics (ASEN6060) - Spring 2024
% Instructor:   Prof. Bosanac
% Student:      Giovanni Fereoli
% Student ID:   111040314

%% Data
clear; clc; close;

% Data
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 
GM_sun = 132712440041.93938;
G_tilde = 6.67408*1e-20; % km^3/kg*s^2

% Computations mass ratios
mu_se = GM_earth / (GM_sun + GM_earth); % Mass Ratio
mu_em = GM_moon / (GM_earth + GM_moon);

%% Problem 1a-1 & 1c
clc; close;

% Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
xi0 = -0.00005;
eta0 = 0;

% Variations and guesses
[initial_var, initial_guess, xL2, P, A] = ...
    init_cond_periodic_planarL2(xi0, eta0, mu_em);

% Correction first guess
tol = 1e-10;
[x0_corr, P_corr, ~] = general_corrector(initial_guess',...
    P, mu_em, tol);

% Continuation scheme 
tol = 1e-12;
step = 0.5 * 1e-2;  
iters = 45; 
[X0_vec, P_vec] = pos_pseudoarc_continuation(x0_corr, P_corr,...  
      step, iters, mu_em, tol, 1, true); 

% Integrate and plot 
disc_plot = 20; 
indices = round(linspace(1, iters, disc_plot));
colormap_vals = linspace(0, 1, disc_plot);
colormap_custom = [colormap_vals.', zeros(disc_plot, 1),...
    1 - colormap_vals.'];
gca1 = figure(1);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
hold on;
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.3 1]);
for i = indices
    %Final integration
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 P_vec(i)],[X0_vec(:,i); reshape(eye(6),[],1)], options);
    % Plot
    norm_x0 = (X0_vec(3, i) - min(X0_vec(3, :))) /...
    (max(X0_vec(3, :)) - min(X0_vec(3, :)));
    color_val = norm_x0;
    plot3(xx(:,1), xx(:,2), xx(:,3), 'Color',...
        colormap_custom(indices == i,:), 'LineWidth', 1.5);
end
quiver3(X0_vec(1,end) + 0.01, X0_vec(2,end), X0_vec(3,end),...
    X0_vec(4,end) / 3, X0_vec(5,end) / 3, X0_vec(6,end) / 3,...
    'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
c = colorbar('eastoutside');
c.Label.String = 'Continuated $x_0$';
c.Label.Interpreter = 'latex';
colormap(c, colormap_custom);
c.TickLabelInterpreter = 'latex'; 
c.Ticks = linspace(0, 1, disc_plot);
c.TickLabels = cellstr(num2str(X0_vec(1, indices).', 3));
grid on;
real axis; 
legend('$L_2$', 'Moon', 'Interpreter', 'latex', 'Location', 'northeast');
view(2);

% Export graphics
exportgraphics(gca1, 'Prob1a1.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Compute stability indeces and plot
stab_vec = stability_index1(X0_vec, P_vec, mu_em);
gca2 = figure(2);
subplot(1, 2, 1);
plot(P_vec(1:10:end), stab_vec(1:10:end, 1), 'r.', 'MarkerSize', 5);
hold on; 
plot(P_vec, 2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
plot(P_vec, -2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
xlabel('$P$ [-]', 'interpreter','latex');
ylabel('$S_1$ [-]', 'interpreter','latex');
xlim([P_vec(1), P_vec(end)]);
real axis;
grid on;
subplot(1, 2, 2);
plot(P_vec(1:10:end), stab_vec(1:10:end, 2), 'b.', 'MarkerSize', 5);
hold on;
plot(P_vec, 2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
plot(P_vec, -2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
xlabel('$P$ [-]', 'interpreter','latex');
ylabel('$S_2$ [-]', 'interpreter','latex');
xlim([P_vec(1), P_vec(end)]);
real axis;
grid on;

% Export graphics
exportgraphics(gca2, 'Prob1a1c.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Exercise 1a-2  & 1c
clc; close;

%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[ ~, xL2, ~, ~, ~] = EquilibriumPoints(mu_em);
tol = 1e-10;
x0guess = [1.180462, 0, -0.0209998, 0, -0.158363, 0]';
P0guess = 3.411921;
[x0_corr, P_corr, ~] = general_corrector(x0guess,...
    P0guess, mu_em, tol);

% Continuation scheme
tol = 1e-10;  
step = 5 * 1e-3;  
iters = 500; 
[X0_vec, P_vec] = neg_pseudoarc_continuation(x0_corr, P_corr,...  
      step, iters, mu_em, tol, 1, true); 

% Integrate and plot 
disc_plot = 30; 
indices = round(linspace(1, iters, disc_plot));
colormap_vals = linspace(0, 1, disc_plot);
colormap_custom = [colormap_vals.', zeros(disc_plot, 1),...
    1 - colormap_vals.'];
gca1 = figure(1);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
hold on;
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.3 1]);
for i = indices
    %Final integration
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 P_vec(i)],[X0_vec(:,i); reshape(eye(6),[],1)], options);
    % Plot
    norm_x0 = (X0_vec(3, i) - min(X0_vec(3, :))) /...
    (max(X0_vec(3, :)) - min(X0_vec(3, :)));
    color_val = norm_x0;
    plot3(xx(:,1), xx(:,2), xx(:,3), 'Color',...
        colormap_custom(indices == i,:), 'LineWidth', 1.5);
end
quiver3(X0_vec(1,end) - 0.03, X0_vec(2,end), X0_vec(3,end),...
    X0_vec(4,end) * 3 , X0_vec(5,end) * 3, X0_vec(6,end) * 3,...
    'k', 'LineWidth', 2);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
c = colorbar('eastoutside');
c.Label.String = 'Continuated $z_0$';
c.Label.Interpreter = 'latex';
colormap(c, colormap_custom);
c.TickLabelInterpreter = 'latex'; 
c.Ticks = linspace(0, 1, disc_plot);
c.TickLabels = cellstr(num2str(X0_vec(3, indices).', 3));
grid on;
real axis; 
legend('$L_2$', 'Moon', 'Interpreter', 'latex', 'Location', 'northeast');

% Export graphics
exportgraphics(gca1, 'Prob1a2.pdf', 'ContentType','image',...
    'Resolution', 1000);

Compute stability indeces and plot
stab_vec = stability_index2(X0_vec, P_vec, mu_em);
gca2 = figure(2);
subplot(1, 2, 1);
plot(P_vec(1:10:end), stab_vec(1:10:end, 1), 'r.', 'MarkerSize', 5);
hold on; 
plot(P_vec, 2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
plot(P_vec, -2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
xlabel('$P$ [-]', 'interpreter','latex');
ylabel('$S_1$ [-]', 'interpreter','latex');
xlim([P_vec(end), P_vec(1)]);
real axis;
grid on;
subplot(1, 2, 2);
plot(P_vec(1:10:end), stab_vec(1:10:end, 2), 'b.', 'MarkerSize', 5);
hold on;
plot(P_vec, 2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
plot(P_vec, -2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
xlabel('$P$ [-]', 'interpreter','latex');
ylabel('$S_2$ [-]', 'interpreter','latex');
xlim([P_vec(end), P_vec(1)]);
real axis;
grid on;

% Export graphics
exportgraphics(gca2, 'Prob1a2c.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Exercise 1a-3  & 1c
clc; close;

% Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[ ~, xL2, ~, ~, ~] = EquilibriumPoints(mu_em);
tol = 1e-3; 
x0guess = [1.0301513, 0, 0, 0, 0.7030025, 0.1552945]';
P0guess = 4.312367;  
[x0_corr, P_corr, ~] = general_corrector(x0guess,...
    P0guess, mu_em, tol);

% Continuation scheme
tol = 1e-10;  
step = 1e-3;  
iters = 1600; 
[X0_vec, P_vec] =  pos_pseudoarc_continuation2(x0_corr, P_corr,...  
      step, iters, mu_em, tol, 1, true); 

% Integrate and plot 
disc_plot = 10; 
indices = round(linspace(1, iters, disc_plot));
colormap_vals = linspace(0, 1, disc_plot);
colormap_custom = [colormap_vals.', zeros(disc_plot, 1),...
    1 - colormap_vals.'];
gca1 = figure(1);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
hold on;
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.3 1]);
for i = indices
    %Final integration
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 P_vec(i)],[X0_vec(:,i); reshape(eye(6),[],1)], options);
    % Plot
    norm_x0 = (X0_vec(3, i) - min(X0_vec(3, :))) /...
    (max(X0_vec(3, :)) - min(X0_vec(3, :)));
    color_val = norm_x0;
    plot3(xx(:,1), xx(:,2), xx(:,3), 'Color',...
        colormap_custom(indices == i,:), 'LineWidth', 1.5);
end
quiver3(X0_vec(1,end) + 0.01, X0_vec(2,end), X0_vec(3,end),...
    X0_vec(4,end) / 2, X0_vec(5,end) / 2, X0_vec(6,end) / 2,...
    'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
c = colorbar('eastoutside');
c.Label.String = 'Continuated $P$';
c.Label.Interpreter = 'latex';
colormap(c, colormap_custom);
c.TickLabelInterpreter = 'latex'; 
c.Ticks = linspace(0, 1, disc_plot);
c.TickLabels = cellstr(num2str(P_vec(indices), 3));
grid on;
real axis; 
legend('$L_2$', 'Moon', 'Interpreter', 'latex', 'Location', 'northeast');

% Export graphics
exportgraphics(gca1, 'Prob1a3.pdf', 'ContentType','image',...
    'Resolution', 1000);

% % Compute stability indeces and plot
% stab_vec = stability_index3(X0_vec, P_vec, mu_em);
% gca2 = figure(2);
% subplot(1, 2, 1);
% plot(P_vec(1:50:end), stab_vec(1:50:end, 1), 'r.', 'MarkerSize', 5);
% hold on; 
% plot(P_vec, 2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
% plot(P_vec, -2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
% xlabel('$P$ [-]', 'interpreter','latex');
% ylabel('$S_1$ [-]', 'interpreter','latex');
% xlim([min(P_vec), max(P_vec)]);
% real axis;
% grid on;
% subplot(1, 2, 2);
% plot(P_vec(1:50:end), stab_vec(1:50:end, 2), 'b.', 'MarkerSize', 5);
% hold on;
% plot(P_vec, 2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
% plot(P_vec, -2*ones(length(P_vec)), 'k--', 'LineWidth', 1);
% xlabel('$P$ [-]', 'interpreter','latex');
% ylabel('$S_2$ [-]', 'interpreter','latex');
% xlim([min(P_vec), max(P_vec)]);
% real axis;
% grid on;
% 
% % Export graphics
% exportgraphics(gca2, 'Prob1a3c.pdf', 'ContentType','image',...
%     'Resolution', 1000);

%% Exercise 1b
clc; close;

%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
tol = 1e-8;
x0guess = [1.180462, 0, -0.0209998, 0, -0.158363, 0]';
P0guess = 3.411921;

% Correction initial guess
[x0_corr, P_corr, ~] = general_corrector(x0guess, P0guess, mu_em, tol);

% Integration with STM
[~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), [0 P_corr], ...
    [x0_corr; reshape(eye(6),[],1)], options);

% Generate Monodromy Matrix
M = (reshape(xx(end, 7:42), 6, 6))';

% Analize monodromy matrix
[eig_vec, eig_val] = eig(M);
disp(eig_vec);
disp(eig_val);

% Accuracy Monodormy Matrix
acc_M = abs(det(M)-1);
disp(acc_M);

%% Exercise 2a/b
clc; close;

%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
options_Moon = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events',@event_Moon);
options_Earth = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events',@event_Earth);
l_star_em = 384400;
d = 15 / l_star_em;
disc = 40;
x0 = [0.826430558525249, 0, 0, 0, 0.095445382297993, 0]';
P0 = 2.720098987628122;

% Generate initial conditions manifolds
[X0_man_stab_pos, X0_man_unstab_pos] = ic_manifolds_pos(x0, P0, d, disc,...
    mu_em);
[X0_man_stab_neg, X0_man_unstab_neg] = ic_manifolds_neg(x0, P0, d, disc,...
    mu_em);

% Integration with STM
[~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), [0 P0], ...
    [x0; reshape(eye(6),[],1)], options);
[~, xx_unstab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 10 * P0], ...
    [X0_man_unstab_pos(:, 1); reshape(eye(6),[],1)], options_Moon);
[~, xx_stab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 -10 * P0], ...
    [X0_man_stab_pos(:, 1); reshape(eye(6),[],1)], options_Moon);
[~, xx_unstab_neg] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 10 * P0], ...
    [X0_man_unstab_neg(:, 1); reshape(eye(6),[],1)], options_Earth);
[~, xx_stab_neg] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 -10 * P0], ...
    [X0_man_stab_neg(:, 1); reshape(eye(6),[],1)], options_Earth);

%Plot trajectory
gca = figure(1);
plot3(xx_stab_pos(:,1), xx_stab_pos(:,2), xx_stab_pos(:,3), 'b',...
    'LineWidth', 0.5);
hold on;
plot3(xx_unstab_pos(:,1), xx_unstab_pos(:,2), xx_unstab_pos(:,3), 'r', ...
    'LineWidth', 0.5);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.3 1]);
plot3(-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.4 0.2 1]);
plot3(xx_stab_neg(:,1), xx_stab_neg(:,2),xx_stab_neg(:,3), 'b',...
    'LineWidth', 0.5);
plot3(xx_unstab_neg(:,1), xx_unstab_neg(:,2),xx_unstab_neg(:,3), 'r', ...
    'LineWidth', 0.5);
for i=2:disc
    % Integration with STM - positive manifolds
    [~, xx_unstab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 10 * P0], ...
        [X0_man_unstab_pos(:, i); reshape(eye(6),[],1)], options_Moon);
    [~, xx_stab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 -10 * P0], ...
        [X0_man_stab_pos(:, i); reshape(eye(6),[],1)], options_Moon);

    % Plot
    plot3(xx_stab_pos(:,1), xx_stab_pos(:,2), xx_stab_pos(:,3), 'b', ...
        'LineWidth', 0.5);
    plot3(xx_unstab_pos(:,1), xx_unstab_pos(:,2), xx_unstab_pos(:,3), 'r',...
        'LineWidth', 0.5);

    % Integration with STM - negative manifolds
    [~, xx_unstab_neg] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 10 * P0], ...
        [X0_man_unstab_neg(:, i); reshape(eye(6),[],1)], options_Earth);
    [~, xx_stab_neg] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), ...
        [0 -10 * P0], ...
        [X0_man_stab_neg(:, i); reshape(eye(6),[],1)], options_Earth);

    % Plot
    plot3(xx_stab_neg(:,1), xx_stab_neg(:,2), xx_stab_neg(:,3), 'b',...
        'LineWidth', 0.5);
    plot3(xx_unstab_neg(:,1), xx_unstab_neg(:,2), xx_unstab_neg(:,3), 'r',...
        'LineWidth', 0.5);
end
plot3(xx(:,1), xx(:,2),xx(:,3), 'k', 'LineWidth', 3);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
view(2);
legend('Stable Manifold', 'Unstable Manifold', 'Moon', 'Earth',...
    'interpreter', 'latex', 'Location','northeast');

% Export graphics
exportgraphics(gca, 'Prob2ab.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Exercise 2c
clc; close;

%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[xL1, ~, ~, ~, ~] = EquilibriumPoints(mu_em);
l_star_em = 384400;
d = 0.1 / l_star_em; 
disc = 50;
x0 = [0.826430558525249, 0, 0, 0, 0.095445382297993, 0]';
P0 = 2.720098987628122;
T_man = 6;

% Generate initial conditions manifolds
[X0_man_stab_pos, X0_man_unstab_pos] = ic_manifolds_pos(x0, P0, d, disc,...
    mu_em);

% Integration with STM
[~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), [0 P0], ...
    [x0; reshape(eye(6),[],1)], options);
[~, xx_unstab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), ...
    [0 T_man], [X0_man_unstab_pos(:, 1); reshape(eye(6),[],1)], options);
[~, xx_stab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), ...
    [0 -T_man], [X0_man_stab_pos(:, 1); reshape(eye(6),[],1)], options);

%Plot trajectory
gca = figure(1);
plot3(xx_stab_pos(:,1), xx_stab_pos(:,2),xx_stab_pos(:,3),'b',...
    'LineWidth', 0.5);
hold on;
plot3(xx_unstab_pos(:,1), xx_unstab_pos(:,2),xx_unstab_pos(:,3),'r',...
    'LineWidth', 0.5);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.3 1]);
for i=2:disc
    % Integration with STM - positive manifolds
    [~, xx_unstab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 T_man], [X0_man_unstab_pos(:, i); reshape(eye(6),[],1)], options);
    [~, xx_stab_pos] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 -T_man], [X0_man_stab_pos(:, i); reshape(eye(6),[],1)], options);

    % Plot
    plot3(xx_stab_pos(:,1), xx_stab_pos(:,2),xx_stab_pos(:,3),'b',...
        'LineWidth', 0.5);
    plot3(xx_unstab_pos(:,1), xx_unstab_pos(:,2),xx_unstab_pos(:,3),'r',...
        'LineWidth', 0.5);
end
plot3(xx(:,1), xx(:,2), xx(:,3),'k','LineWidth', 3);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
view(2);
legend('Stable Manifold', 'Unstable Manifold',...
    '$L_1$', 'Moon', 'interpreter', 'latex', 'Location',...
    'northeast');

% Export graphics
exportgraphics(gca, 'Prob2c.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Functions

% Find equilibrium points
function [xL1, xL2, xL3, xL4, xL5] = EquilibriumPoints(mu)
% Collinear points
% Position primaries along x
xxP1 = -mu;
xxP2 = 1-mu;

% Gradient of U* in x
f = @(x) x-(1-mu)*(x+mu)/(abs(x+mu))^3-mu*(x+mu-1)/(abs(x+mu-1))^3; 

%Inital guesses
z = (mu/3)^(1/3);
xxL10 = xxP2 - (z-(1/3)*z^2-(1/9)*z^3+(58/81)*z^4);
xxL20 = xxP2 + (z+(1/3)*z^2-(1/9)*z^3+(50/81)*z^4);
xxL30 = xxP1 - (1-(7/12)*mu-(1127/20736)*mu^3-(7889/248832)*mu^4);

%Zeros computation
options = optimoptions('fsolve', 'Display', 'none', 'TolFun', 1e-15);
xL1 = [fsolve(f, xxL10, options), 0, 0];
xL2 = [fsolve(f, xxL20, options), 0, 0]; 
xL3 = [fsolve(f, xxL30, options), 0, 0];

% Triangular points
xL4 = [0.5 - mu, sqrt(3) / 2, 0];
xL5 = [0.5 - mu, - sqrt(3) / 2, 0];

end

% C3RBP ODEs: state and STM
function dxMdt = STMexact_CRTBP(~,xM,mu)
%Initialize ODE
dxMdt = zeros(42,1);

% Unpack state
x = xM(1);
y = xM(2);
z = xM(3);
xdot = xM(4);
ydot = xM(5);
zdot = xM(6);

% Unpack STM
M = (reshape(xM(7:42),6,6))';   %From state to STM

% CR3TBP miscellaneous
% P3 distance from primaries
r1_norm = sqrt((x+mu)^2+y^2+z^2);
r2_norm = sqrt((x+mu-1)^2+y^2+z^2);

% Variational equations
df4dx = 1-(1-mu)/r1_norm^3+3*(1-mu)*(x+mu)^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*(x+mu-1)^2/r2_norm^5;
df4dy = 3*(1-mu)*(x+mu)*y/r1_norm^5+3*mu*(x+mu-1)*y/r2_norm^5;
df4dz = 3*(1-mu)*(x+mu)*z/r1_norm^5+3*mu*(x+mu-1)*z/r2_norm^5;
df5dy = 1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;
df5dz = 3*(1-mu)*y*z/r1_norm^5+3*mu*y*z/r2_norm^5;
df6dz = -(1-mu)/r1_norm^3+3*(1-mu)*z^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*z^2/r2_norm^5;

% Jacobian 
A = [0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 1;...
    df4dx, df4dy, df4dz, 0, 2, 0;...
    df4dy, df5dy, df5dz, -2, 0, 0;...
    df4dz, df5dz, df6dz, 0, 0, 0];

% CR3TBP dynamics
% CR3BP dynamics: position and velocity
dxMdt(1:3) = [xdot; ydot; zdot];
dxMdt(4:6) = [2*ydot+x-(1-mu)*(x+mu)/r1_norm^3-mu*(x+mu-1)/r2_norm^3;...
    -2*xdot+y-(1-mu)*y/r1_norm^3-mu*y/r2_norm^3;...
    -(1-mu)*z/r1_norm^3-mu*z/r2_norm^3];

% CR3BP dynamics: STM
dMdt = A*M;
dxMdt(7:12) = dMdt(1,1:6)';
dxMdt(13:18) = dMdt(2,1:6)';
dxMdt(19:24) = dMdt(3,1:6)';
dxMdt(25:30) = dMdt(4,1:6)';
dxMdt(31:36) = dMdt(5,1:6)';
dxMdt(37:42) = dMdt(6,1:6)';

end

% General initial guess correction algorithm
function [X0_corr, P_corr, normF_iter] = general_corrector(X0, P, mu, tol)
% Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
F = ones(6,1);
V = [X0; P];
normF_iter = zeros(100, 1);
iter = 1;

while norm(F) > tol
    %Integration state and STM
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)],...
    [V(1:6); reshape(eye(6), [], 1)], options);

    %Reshape final STM
    STMf=(reshape(xx(end, 7:42), 6, 6))';

    %Jacobian of F
    xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
    DF = [STMf - eye(6), xMdotf(1:6)];

    % Update F
    F = xx(end,1:6)- xx(1,1:6);
    normF_iter(iter) = norm(F);

    % Check if norm of F exceeds tolerance
    if norm(F) > 1e3
        fprintf('Error: Norm of F diverging. \n');
        break;  % Exit the loop if norm(F) exceeds tolerance
    end

    % Correction
    dV = - lsqminnorm(DF, F', tol);
    V = V + dV;
    iter = iter + 1;
end

% Final correction
X0_corr = V(1:6);
P_corr = V(7);
normF_iter = normF_iter(1:iter-1);

% Final print
disp('Correction terminated successfully.');
fprintf('         Current function value: %.14f\n', normF_iter(end));
fprintf('         Iterations: %d\n', iter-1);
end

% Pseudo-Arclength continuation positive
function [X0_vec, P_vec] = pos_pseudoarc_continuation(X0_in, P_in, ...
    step, iters, mu, tol, comp_cont, print)
% Initialization
if nargin < 8
    print = false;
end
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
X0_vec = zeros(6, iters);
P_vec = zeros(iters, 1);
iter_cont = 1;
V = [X0_in; P_in];
X0_vec(:, 1) = V(1:6);  % OSS: here, change only P and not V a-priori.
P_vec(1) = V(7);
mat_mod = [eye(4), zeros(4,2); zeros(1,5), 1];

while iter_cont < iters
    % Initilization first loop
    H = ones(7,1); 
    iter_corr = 1;
    flag = 0;

    % Pseudo-arclength initial guess with checks for x0-continuation
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)], ...
        [V(1:6); reshape(eye(6), [], 1)], options);
    STMf = (reshape(xx(end, 7:42), 6, 6))';
    xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
    DF = [STMf([1:4, 6], :) - mat_mod, xMdotf([1:4, 6]);...
        0, 1, 0, 0, 0, 0, 0];
    nullDF = null(DF, tol);
    nullDF = nullDF(:, 1);   % Check columns 
    if nullDF(comp_cont) < 0
        nullDF = - nullDF;
    end
    Vd = V;
    V = Vd + step * nullDF; % Initial guess of the next orbit

    while norm(H) > tol
        %Integration state and STM
        [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)],...
            [V(1:6); reshape(eye(6), [], 1)], options);
    
        %Reshape final STM
        STMf = (reshape(xx(end, 7:42), 6, 6))';
    
        %Jacobian of H
        xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
        DH = [STMf([1:4, 6], :) - mat_mod, xMdotf([1:4, 6]);...
            0, 1, 0, 0, 0, 0, 0;...
            nullDF'];
    
        % Update H
        H(1:5) = xx(end, [1:4, 6])- xx(1, [1:4, 6]);
        H(6) = xx(1, 2);
        H(7) = dot(V-Vd, nullDF) - step;

        % Check if norm of H exceeds tolerance
        if norm(H) > 1e3
            fprintf('Error: Norm of H diverging.');
            flag = 1;
            break;  % Exit the loop if norm(H) exceeds tolerance
        end
    
        % Correction
        dV = - lsqminnorm(DH, H, tol);
        V = V + dV;
        iter_corr = iter_corr + 1;
    end

    % Check 
    if flag == 1
        break
    end

    % Final print
    if print == true
        disp('Correction terminated successfully.');
        fprintf('         Current function value: %.14f\n', norm(H));
        fprintf('         Iterations: %d\n', iter_corr);
    end

    % Final correction and updates
    iter_cont = iter_cont + 1;
    X0_vec(:, iter_cont) = V(1:6);
    P_vec(iter_cont) = V(7);
end
end

% Pseudo-Arclength continuation positive
function [X0_vec, P_vec] = neg_pseudoarc_continuation(X0_in, P_in, ...
    step, iters, mu, tol, comp_cont, print)
% Initialization
if nargin < 8
    print = false;
end
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
X0_vec = zeros(6, iters);
P_vec = zeros(iters, 1);
iter_cont = 1;
V = [X0_in; P_in];
X0_vec(:, 1) = V(1:6);  % OSS: here, change only P and not V a-priori.
P_vec(1) = V(7);
mat_mod = [eye(4), zeros(4,2); zeros(1,5), 1];

while iter_cont < iters
    % Initilization first loop
    H = ones(7,1); 
    iter_corr = 1;
    flag = 0;

    % Pseudo-arclength initial guess with checks for x0-continuation
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)], ...
        [V(1:6); reshape(eye(6), [], 1)], options);
    STMf = (reshape(xx(end, 7:42), 6, 6))';
    xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
    DF = [STMf([1:4, 6], :) - mat_mod, xMdotf([1:4, 6]);...
        0, 1, 0, 0, 0, 0, 0];
    nullDF = null(DF, tol);
    nullDF = nullDF(:, 1);   % Check columns 
    if nullDF(comp_cont) > 0
        nullDF = - nullDF;
    end
    Vd = V;
    V = Vd + step * nullDF; % Initial guess of the next orbit

    while norm(H) > tol
        %Integration state and STM
        [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)],...
            [V(1:6); reshape(eye(6), [], 1)], options);
    
        %Reshape final STM
        STMf = (reshape(xx(end, 7:42), 6, 6))';
    
        %Jacobian of H
        xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
        DH = [STMf([1:4, 6], :) - mat_mod, xMdotf([1:4, 6]);...
            0, 1, 0, 0, 0, 0, 0;...
            nullDF'];
    
        % Update H
        H(1:5) = xx(end, [1:4, 6])- xx(1, [1:4, 6]);
        H(6) = xx(1, 2);
        H(7) = dot(V-Vd, nullDF) - step;

        % Check if norm of H exceeds tolerance
        if norm(H) > 1e3
            fprintf('Error: Norm of H diverging.');
            flag = 1;
            break;  % Exit the loop if norm(H) exceeds tolerance
        end
    
        % Correction
        dV = - lsqminnorm(DH, H, tol);
        V = V + dV;
        iter_corr = iter_corr + 1;
    end

    % Check 
    if flag == 1
        break
    end

    % Final print
    if print == true
        disp('Correction terminated successfully.');
        fprintf('         Current function value: %.14f\n', norm(H));
        fprintf('         Iterations: %d\n', iter_corr);
    end

    % Final correction and updates
    iter_cont = iter_cont + 1;
    X0_vec(:, iter_cont) = V(1:6);
    P_vec(iter_cont) = V(7);
end
end


% Compute analytically initial guess for Periodic Orbits about L2
function [initial_var, initial_guess, xL2, P, A] =...
    init_cond_periodic_planarL2(xi0, eta0, mu)

% Compute equilibrium point, Jacobian and eigenvalues
[~, xL2, ~, ~, ~] = EquilibriumPoints(mu);
A = Jacobian_CRTBP(xL2, mu);
[eig_L2_inplane, ~] = Eigenvalues_CRTBP_Lagrange(xL2, mu);

% Miscellaneous
lam3sq = (eig_L2_inplane(3))^2;
Uxx_eq = A(4,1);
P = 2*pi / abs(imag(eig_L2_inplane(3)));

% Compute initial variations
xidot0 = 2 * lam3sq * eta0 / (lam3sq - Uxx_eq);
etadot0 = 0.5 * (lam3sq - Uxx_eq) * xi0;
initial_var = [xi0, eta0, 0, xidot0, etadot0, 0];

% Compute non-linear initial guess
initial_guess = [xL2(1), xL2(2), xL2(3), 0, 0, 0] + initial_var;
end

% Compute Eigenvalues in-plane and out-of-plane, Lagrange points
function [eig_inplane, eig_outplane] = Eigenvalues_CRTBP_Lagrange(X, mu)
% HP: It's for equilibrium points, being in the plane they have
% df1dz=df2dz=0. 

%Initialize
x = X(1);
y = X(2);
z = X(3);
r1_norm = sqrt((x+mu)^2+y^2+z^2);
r2_norm = sqrt((x+mu-1)^2+y^2+z^2);

% Variational equations
df1dx = 1-(1-mu)/r1_norm^3+3*(1-mu)*(x+mu)^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*(x+mu-1)^2/r2_norm^5;
df1dy = 3*(1-mu)*(x+mu)*y/r1_norm^5+3*mu*(x+mu-1)*y/r2_norm^5;
% df1dz = 3*(1-mu)*(x+mu)*z/r1_norm^5+3*mu*(x+mu-1)*z/r2_norm^5;
df2dy = 1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;
% df2dz = 3*(1-mu)*y*z/r1_norm^5+3*mu*y*z/r2_norm^5;
df3dz = -(1-mu)/r1_norm^3+3*(1-mu)*z^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*z^2/r2_norm^5;

% Jacobian
A_inplane = [0, 0, 1, 0;...
    0, 0, 0, 1;...
    df1dx, df1dy, 0, 2;...
    df1dy, df2dy, -2, 0];

A_outplane = [0, 1;...
    df3dz, 0];

% Eigenvalues
eig_inplane = eig(A_inplane);
eig_outplane = eig(A_outplane);

end

% Compute Jacobian of the CR3BP
function A = Jacobian_CRTBP(X, mu)
%Initialize
x = X(1);
y = X(2);
z = X(3);
r1_norm = sqrt((x+mu)^2+y^2+z^2);
r2_norm = sqrt((x+mu-1)^2+y^2+z^2);

% Variational equations
df1dx = 1-(1-mu)/r1_norm^3+3*(1-mu)*(x+mu)^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*(x+mu-1)^2/r2_norm^5;
df1dy = 3*(1-mu)*(x+mu)*y/r1_norm^5+3*mu*(x+mu-1)*y/r2_norm^5;
df1dz = 3*(1-mu)*(x+mu)*z/r1_norm^5+3*mu*(x+mu-1)*z/r2_norm^5;
df2dy = 1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;
df2dz = 3*(1-mu)*y*z/r1_norm^5+3*mu*y*z/r2_norm^5;
df3dz = -(1-mu)/r1_norm^3+3*(1-mu)*z^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*z^2/r2_norm^5;

% Jacobian
A= [0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 1;...
    df1dx, df1dy, df1dz, 0, 2, 0;...
    df1dy, df2dy, df2dz, -2, 0, 0;...
    df1dz, df2dz, df3dz, 0, 0, 0];

end

% Compute stability indeces - Pt1
function stab_vec = stability_index1(x0_vec, P_vec, mu)

% Initialization
disc = length(x0_vec);
stab_vec = zeros(disc, 2);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);

% For loop
for i=1:disc
    % Integration with STM
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 P_vec(i)], ...
        [x0_vec(:, i); reshape(eye(6),[],1)], options);

    % Generate Monodromy Matrix
    M = (reshape(xx(end, 7:42), 6, 6))';

    % Analize monodromy matrix
    [~, eig_val] = eig(M);

    % Find the trivial pair and remove them
    [~, trivial_indices] = sort(abs(diag(eig_val) - 1));
    trivial_indices = trivial_indices(1:2);
    eig_val = eig_val(~ismember(1:numel(diag(eig_val)), trivial_indices),...
        ~ismember(1:numel(diag(eig_val)), trivial_indices));

    % Stability indeces
    % HP: stability index should be Real!
    stab_vec(i, 1) = eig_val(1,1) + eig_val(2,2); 
    stab_vec(i, 2) = eig_val(3,3) + eig_val(4,4); 
end
end

% Compute stability indeces - Pt2
function stab_vec = stability_index2(x0_vec, P_vec, mu)

% Initialization
disc = length(x0_vec);
stab_vec = zeros(disc, 2);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);

% For loop
for i=1:disc
    % Integration with STM
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 P_vec(i)], ...
        [x0_vec(:, i); reshape(eye(6),[],1)], options);

    % Generate Monodromy Matrix
    M = (reshape(xx(end, 7:42), 6, 6))';

    % Analize monodromy matrix
    [~, eig_val] = eig(M);

    % Find the trivial pair and remove them
    [~, trivial_indices] = sort(abs(diag(eig_val) - 1));
    trivial_indices = trivial_indices(1:2);
    eig_val = eig_val(~ismember(1:numel(diag(eig_val)), trivial_indices),...
        ~ismember(1:numel(diag(eig_val)), trivial_indices));

    % Stability indeces
    % HP: stability index should be Real!
    stab_vec(i, 1) = eig_val(1,1) + eig_val(2,2); 
    stab_vec(i, 2) = eig_val(3,3) + eig_val(4,4); 
    if ~isreal(stab_vec(i, 1)) || ~isreal(stab_vec(i, 2)) % Check 1
         stab_vec(i, 1) = eig_val(1,1) + eig_val(3,3);  
         stab_vec(i, 2) = eig_val(2,2) + eig_val(4,4); 
    end
    if ~isreal(stab_vec(i, 1)) || ~isreal(stab_vec(i, 2)) % Check 2
         stab_vec(i, 1) = eig_val(1,1) + eig_val(4,4);  
         stab_vec(i, 2) = eig_val(2,2) + eig_val(3,3); 
    end

    % Check for swapping
    if i>1 && abs(stab_vec(i,2) - stab_vec(i-1, 2)) > 0.1 
         iter_swap = stab_vec(i,2);
         stab_vec(i,2) = stab_vec(i,1);
         stab_vec(i,1) = iter_swap;
    end
end
end

% Compute stability indeces - Pt3
function stab_vec = stability_index3(x0_vec, P_vec, mu)

% Initialization
disc = length(x0_vec);
stab_vec = zeros(disc, 2);
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);

% For loop
for i=1:disc
    % Integration with STM
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 P_vec(i)], ...
        [x0_vec(:, i); reshape(eye(6),[],1)], options);

    % Generate Monodromy Matrix
    M = (reshape(xx(end, 7:42), 6, 6))';

    % Analize monodromy matrix
    [~, eig_val] = eig(M);

    % Find the trivial pair and remove them
    [~, trivial_indices] = sort(abs(diag(eig_val) - 1));
    trivial_indices = trivial_indices(1:2);
    eig_val = eig_val(~ismember(1:numel(diag(eig_val)), trivial_indices),...
        ~ismember(1:numel(diag(eig_val)), trivial_indices));

    % Stability indeces 
    % HP: stability index should be Real!
    if eig_val(1,1) == 1 / eig_val(2,2)
        stab_vec(i, 1) = eig_val(1,1) + eig_val(2,2); 
        stab_vec(i, 2) = eig_val(3,3) + eig_val(4,4); 
    elseif eig_val(1,1) == 1 / eig_val(3,3)
         stab_vec(i, 1) = eig_val(1,1) + eig_val(3,3);  
         stab_vec(i, 2) = eig_val(2,2) + eig_val(4,4); 
    else
         stab_vec(i, 1) = eig_val(1,1) + eig_val(4,4);  
         stab_vec(i, 2) = eig_val(2,2) + eig_val(3,3); 
    end
end
end
 
% % Compute stability indeces - Pt4
% function stab_vec = stability_index4(x0_vec, P_vec, mu)
% 
% % Initialization
% disc = length(x0_vec);
% stab_vec = zeros(disc, 2);
% options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
% 
% % Function for reciprocal, Real and Complex
% function result = reciprocal(x)
%     if isreal(x)
%         result = 1 / x;
%     else
%         result = conj(x);
%     end
% end
% 
% % For loop
% for i=1:disc
%     % Integration with STM
%     [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 P_vec(i)], ...
%         [x0_vec(:, i); reshape(eye(6),[],1)], options);
% 
%     % Generate Monodromy Matrix
%     M = (reshape(xx(end, 7:42), 6, 6))';
% 
%     % Analize monodromy matrix
%     [~, eig_val] = eig(M);
% 
%     % Find the trivial pair and remove them
%     [~, trivial_indices] = sort(abs(diag(eig_val) - 1));
%     trivial_indices = trivial_indices(1:2);
%     eig_val = eig_val(~ismember(1:numel(diag(eig_val)), trivial_indices),...
%         ~ismember(1:numel(diag(eig_val)), trivial_indices));
% 
%     % Stability indeces 
%     % HP: stability index should be Real!
%     if abs(reciprocal(eig_val(1,1)) - eig_val(2,2)) < 1e-8
%         stab_vec(i, 1) = eig_val(1,1) + eig_val(2,2); 
%         stab_vec(i, 2) = eig_val(3,3) + eig_val(4,4); 
%     elseif abs(reciprocal(eig_val(1,1)) - eig_val(3,3)) < 1e-8
%          stab_vec(i, 1) = eig_val(1,1) + eig_val(3,3);  
%          stab_vec(i, 2) = eig_val(2,2) + eig_val(4,4); 
%     else
%          stab_vec(i, 1) = eig_val(1,1) + eig_val(4,4);  
%          stab_vec(i, 2) = eig_val(2,2) + eig_val(3,3); 
%     end
% end
% end

% Pseudo-Arclength continuation positive - Pt2
function [X0_vec, P_vec] = pos_pseudoarc_continuation2(X0_in, P_in, ...
    step, iters, mu, tol, comp_cont, print)
% Initialization
if nargin < 8
    print = false;
end
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
X0_vec = zeros(6, iters);
P_vec = zeros(iters, 1);
iter_cont = 1;
V = [X0_in; P_in];
X0_vec(:, 1) = V(1:6);  % OSS: here, change only P and not V a-priori.
P_vec(1) = V(7);
mat_mod = [eye(4), zeros(4,2); zeros(1,4), 1, 0];

while iter_cont < iters
    % Initilization first loop
    H = ones(7,1); 
    iter_corr = 1;
    flag = 0;

    % Pseudo-arclength initial guess with checks for x0-continuation
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)], ...
        [V(1:6); reshape(eye(6), [], 1)], options);
    STMf = (reshape(xx(end, 7:42), 6, 6))';
    xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
    DF = [STMf(1:5, :) - mat_mod, xMdotf(1:5);...
        0, 0, 1, 0, 0, 0, 0];
    nullDF = null(DF, tol);
    nullDF = nullDF(:, 1);   % Check columns 
    if nullDF(comp_cont) < 0
        nullDF = - nullDF;
    end
    Vd = V;
    V = Vd + step * nullDF; % Initial guess of the next orbit

    while norm(H) > tol
        %Integration state and STM
        [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)],...
            [V(1:6); reshape(eye(6), [], 1)], options);
    
        %Reshape final STM
        STMf = (reshape(xx(end, 7:42), 6, 6))';
    
        %Jacobian of H
        xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
        DH = [STMf(1:5, :) - mat_mod, xMdotf(1:5);...
            0, 0, 1, 0, 0, 0, 0;...
            nullDF'];
    
        % Update H
        H(1:5) = xx(end, 1:5) - xx(1, 1:5);
        H(6) = xx(1, 3);
        H(7) = dot(V-Vd, nullDF) - step;

        % Check if norm of H exceeds tolerance
        if norm(H) > 1e3
            fprintf('Error: Norm of H diverging.');
            flag = 1;
            break;  % Exit the loop if norm(H) exceeds tolerance
        end
    
        % Correction
        dV = - lsqminnorm(DH, H, tol);
        V = V + dV;
        iter_corr = iter_corr + 1;
    end

    % Check 
    if flag == 1
        break
    end

    % Final print
    if print == true
        disp('Correction terminated successfully.');
        fprintf('         Current function value: %.14f\n', norm(H));
        fprintf('         Iterations: %d\n', iter_corr);
    end

    % Final correction and updates
    iter_cont = iter_cont + 1;
    X0_vec(:, iter_cont) = V(1:6);
    P_vec(iter_cont) = V(7);
end
end

% Positive half-manifold
function [X0_man_stab, X0_man_unstab] = ic_manifolds_pos(x0, P, d, disc, mu)
% Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), ...
    linspace(0, P, disc), [x0; reshape(eye(6),[],1)], options);
X0_man_unstab = zeros(6, disc);
X0_man_stab = zeros(6, disc);

% Generate x0 Monodromy Matrix
M = (reshape(xx(end, 7:42), 6, 6))';

% Analize monodromy matrix
[eig_vec, ~] = eig(M);

% Generate initial directions
v_unstab_manifold = eig_vec(:,1);  % HP: always first and second!
v_stab_manifold = eig_vec(:,2);

% Initial IC's manifolds ICs
X0_man_unstab(:, 1)  = x0 + d * v_unstab_manifold / norm(x0(1:3));
X0_man_stab(:, 1) =  x0 + d * v_stab_manifold / norm(x0(1:3));

% Loop for all manifold' ICs
for i = 2:disc
    % Generate state and STM iteration
    STM_iter = (reshape(xx(i, 7:42), 6, 6))';
    xx_iter = xx(i, 1:6)';

    % Initial IC's manifolds ICs
    X0_man_unstab(:, i)  = xx_iter + d * STM_iter * v_unstab_manifold...
        / norm(STM_iter * v_unstab_manifold);
    X0_man_stab(:, i) =  xx_iter + d * STM_iter * v_stab_manifold...
        / norm(STM_iter * v_stab_manifold);
end
end

% Negative half-manifold
function [X0_man_stab, X0_man_unstab] = ic_manifolds_neg(x0, P, d, disc, mu)
% Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), ...
    linspace(0, P, disc), [x0; reshape(eye(6),[],1)], options);
X0_man_unstab = zeros(6, disc);
X0_man_stab = zeros(6, disc);

% Generate x0 Monodromy Matrix
M = (reshape(xx(end, 7:42), 6, 6))';

% Analize monodromy matrix
[eig_vec, ~] = eig(M);

% Generate initial directions
v_unstab_manifold = - eig_vec(:,1);  % HP: always first and second!
v_stab_manifold = - eig_vec(:,2);

% Initial IC's manifolds ICs
X0_man_unstab(:, 1)  = x0 + d * v_unstab_manifold / norm(x0(1:3));
X0_man_stab(:, 1) =  x0 + d * v_stab_manifold / norm(x0(1:3));

% Loop for all manifold' ICs
for i = 2:disc
    % Generate state and STM iteration
    STM_iter = (reshape(xx(i, 7:42), 6, 6))';
    xx_iter = xx(i, 1:6)';

    % Initial IC's manifolds ICs
    X0_man_unstab(:, i)  = xx_iter + d * STM_iter * v_unstab_manifold...
        / norm(STM_iter * v_unstab_manifold);
    X0_man_stab(:, i) =  xx_iter + d * STM_iter * v_stab_manifold...
        / norm(STM_iter * v_stab_manifold);
end
end

% Event Moon
function [value, isterminal, direction] = event_Moon(~, X, ~)
% Data
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 

% Computations mass ratios
mu_em = GM_moon / (GM_earth + GM_moon);

% Event
value = X(1) - (1-mu_em);   % x = 1-mu
isterminal = 1;
direction = 1; % all directions
end


% Event Earth
function [value, isterminal, direction] = event_Earth(~, X, ~)
% Data
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 

% Computations mass ratios
mu_em = GM_moon / (GM_earth + GM_moon);

%Event
value = X(1) + mu_em;   % x = -mu
isterminal = 1;
direction = 0; % all directions
end







