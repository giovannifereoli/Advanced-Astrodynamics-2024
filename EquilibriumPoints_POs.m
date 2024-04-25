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

%% Problem 1

% Computation Equilibrium Points
[xL1, xL2, xL3, xL4, xL5] = EquilibriumPoints(mu_em);
C_L1 = JacobiConstant([xL1, zeros(1, 3)], mu_em);
C_L2 = JacobiConstant([xL2, zeros(1, 3)], mu_em);
C_L3 = JacobiConstant([xL3, zeros(1, 3)], mu_em);
C_L4 = JacobiConstant([xL4, zeros(1, 3)], mu_em);
C_L5 = JacobiConstant([xL5, zeros(1, 3)], mu_em);

% Effect of mass ratio on the equilibrium points
disc = 1000;
mu_vec = linspace(mu_se, 0.4, disc);
pos_lag_mu = zeros(disc, 5, length(xL1));
P1 = zeros(disc);
P2 = zeros(disc);

for j = 1:disc
    [xL1_step, xL2_step, xL3_step,...
        xL4_step, xL5_step] = EquilibriumPoints(mu_vec(j));
    P1(j) = -mu_vec(j);
    P2(j) = 1-mu_vec(j);
    pos_lag_mu(j, 1, :) = xL1_step;
    pos_lag_mu(j, 2, :) = xL2_step;
    pos_lag_mu(j, 3, :) = xL3_step;
    pos_lag_mu(j, 4, :) = xL4_step;
    pos_lag_mu(j, 5, :) = xL5_step;
end

figure(1);
scatter(pos_lag_mu(:, 1, 1), pos_lag_mu(:, 1, 2), 50, mu_vec, 'filled');
text(pos_lag_mu(round(disc/3), 1, 1), pos_lag_mu(1, 1, 2)+0.15, '$L_1$',...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
hold on;
scatter(pos_lag_mu(:, 2, 1), pos_lag_mu(:, 2, 2), 50, mu_vec, 'filled');
text(pos_lag_mu(round(disc/20), 2, 1), pos_lag_mu(1, 2, 2)+0.15, '$L_2$',...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
scatter(pos_lag_mu(:, 3, 1), pos_lag_mu(:, 3, 2), 50, mu_vec, 'filled');
text(pos_lag_mu(round(disc/2), 3, 1), pos_lag_mu(1, 3, 2)+0.15, '$L_3$',...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
scatter(pos_lag_mu(:, 4, 1), pos_lag_mu(:, 4, 2), 50, mu_vec, 'filled');
text(pos_lag_mu(round(disc/2), 4, 1), pos_lag_mu(1, 4, 2)+0.15, '$L_4$',...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
scatter(pos_lag_mu(:, 5, 1), pos_lag_mu(:, 5, 2), 50, mu_vec, 'filled');
text(pos_lag_mu(round(disc/2), 5, 1), pos_lag_mu(1, 5, 2)+0.15, '$L_5$',...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
colormap(jet);
colorbar;
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
grid on;
ax = gca;
ax_pos = ax.Position;
cbar = colorbar;
title (cbar, '$\mu$ [-]', 'interpreter','latex');

% Export graphics
exportgraphics(gca, 'Prob1C.pdf', 'ContentType','image',...
    'Resolution', 1000);
close all;

% Plot 2
figure(2);
plot(mu_vec, pos_lag_mu(:, 1, 1), 'r', LineWidth=1.5);
hold on;
plot(mu_vec, pos_lag_mu(:, 2, 1), 'b', LineWidth=1.5);
plot(mu_vec, pos_lag_mu(:, 3, 1), 'g', LineWidth=1.5);
plot(mu_vec, pos_lag_mu(:, 4, 1), 'k', LineWidth=1.5);
plot(mu_vec, pos_lag_mu(:, 5, 1), 'c', LineWidth=1.5);
plot(mu_vec, -mu_vec, '--', LineWidth=1.5);
plot(mu_vec, 1-mu_vec, '--', LineWidth=1.5);
xlabel('$\mu$ [-]', 'interpreter','latex');
ylabel('x [-]', 'interpreter','latex');
legend('$x_{L_1}$', '$x_{L_2}$', '$x_{L_3}$', '$x_{L_4}$', '$x_{L_5}$',...
  '$x_{P_1}$', '$x_{P_2}$', 'interpreter','latex', 'Location', 'northeast');
real axis;
grid on;

% Export graphics 2
exportgraphics(gca, 'Prob1Cb.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 2
clc; close;

% Format
format bank;

% Equilibrium points
[xL1_em, xL2_em, xL3_em, xL4_em, xL5_em] = EquilibriumPoints(mu_em);
[xL1_se, xL2_se, xL3_se, xL4_se, xL5_se] = EquilibriumPoints(mu_se);

% Eigenvalues equilibrium points Earth-Moon
[eig_L1_em_inplane, eig_L1_em_outplane] = Eigenvalues_CRTBP_Lagrange(xL1_em, mu_em);
[eig_L2_em_inplane, eig_L2_em_outplane] = Eigenvalues_CRTBP_Lagrange(xL2_em, mu_em);
[eig_L3_em_inplane, eig_L3_em_outplane] = Eigenvalues_CRTBP_Lagrange(xL3_em, mu_em);
[eig_L4_em_inplane, eig_L4_em_outplane] = Eigenvalues_CRTBP_Lagrange(xL4_em, mu_em);
[eig_L5_em_inplane, eig_L5_em_outplane] = Eigenvalues_CRTBP_Lagrange(xL5_em, mu_em);

% Eigenvalues equilibrium points Sun-Earth
[eig_L1_se_inplane, eig_L1_se_outplane] = Eigenvalues_CRTBP_Lagrange(xL1_se, mu_se);
[eig_L2_se_inplane, eig_L2_se_outplane] = Eigenvalues_CRTBP_Lagrange(xL2_se, mu_se);
[eig_L3_se_inplane, eig_L3_se_outplane] = Eigenvalues_CRTBP_Lagrange(xL3_se, mu_se);
[eig_L4_se_inplane, eig_L4_se_outplane] = Eigenvalues_CRTBP_Lagrange(xL4_se, mu_se);
[eig_L5_se_inplane, eig_L5_se_outplane] = Eigenvalues_CRTBP_Lagrange(xL5_se, mu_se);


%% Problem 3

% Initialization
disc = 10000;
mu_vec = linspace(0.0383, 0.039, disc);
real_eigmax_L5 = zeros(disc, 1);
mu_crit = 0;

% Effect of mass ratio on L5 stability
for j=1:disc
   % Location L5
   [~, ~, ~, ~, xL5] = EquilibriumPoints(mu_vec(j));

   % Max real part eigenvalues of Jacobian
   real_eigmax_L5(j) = max(round(real(eig(Jacobian_CRTBP(xL5,...
       mu_vec(j)))), 15));

   % Save critical mass ration
   if real_eigmax_L5(j)>0
       if mu_crit == 0
           mu_crit = mu_vec(j);
       end
   end
end

% Plot
figure(1);
semilogy(mu_vec, real_eigmax_L5, 'r', 'LineWidth', 1);
xlabel('$\mu$ [-]', 'interpreter','latex');
ylabel('$max(Real(\lambda_{1:4}))$ [-]', 'interpreter','latex');
grid on;
ax = gca;

% Export graphics
exportgraphics(gca, 'Prob3.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 4C
clc; close;

% Initialization
disc = 250;
xi0 = -0.00005;
eta0 = 0;

% Variations and guesses
[initial_var, initial_guess, xL1, P, A] = ...
    init_cond_periodic_planarL1(xi0, eta0, mu_em);

% Integrate variation CR3BP and obtain trajectory
tspan = linspace(0, P, disc);
options = odeset('RelTol', 2.22045*1e-14, 'AbsTol', 2.22045*1e-16);
[~, xx1_var] = ode113(@(t,X) CRTBP_lin(t, X, A), tspan, initial_var, options);
xx1_lin = xx1_var;
for j = 1:length(xx1_var)
    xx1_lin(j, 1:3) = xx1_lin(j, 1:3) + xL1;
end

% Integrate CR3BP
[~, xx1] = ode113(@(t,X) CRTBP(t, X, mu_em), tspan, initial_guess, options);
 
% Plot
gca = figure(1);
plot3(xx1_lin(:,1), xx1_lin(:,2), xx1_lin(:,3), 'b', 'LineWidth', 1);
hold on;
plot3(xx1(:,1), xx1(:,2), xx1(:,3), 'r', 'LineWidth', 1);
plot3(xL1(1), xL1(2), xL1(3), 'c.', 'MarkerSize', 10);
plot3(xx1(1,1), xx1(1,2), xx1(1,3), 'k*', 'Markersize', 7);
quiver3(xx1(1,1)-0.0004, xx1(1,2), xx1(1,3), xx1(1,4)/4, xx1(1,5)/4,...
    xx1(1,6)/4, 'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
legend('Linearized CR3BP', 'CR3BP', '$L_1$', '$\mathbf{x}_0$',...
    'interpreter','latex',...
    'Location','southwest', 'FontSize', 10);
grid on;
real axis;
view(2);

% Save plot
exportgraphics(gca, 'Prob4C.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot
norm_diff = zeros(disc,1);
for i = 1:disc
    norm_diff(i) = norm(xx1_lin(i, :) - xx1(i, :));
end
gca = figure(2);
semilogy(tspan, norm_diff, 'k', 'LineWidth', 1);
xlabel('t [-]', 'interpreter','latex');
ylabel('$||\mathbf{x}-\mathbf{x}_{lin}||$ [-]', 'interpreter','latex');
grid on;
real axis;
view(2);

% Save plot
exportgraphics(gca, 'Prob4Cdistance.pdf', 'ContentType','image',...
    'Resolution', 1000);


%% Problem 4D
clc; close;

% L1 eigenvectors
[eigenvec_L1, ~] = Eig_CRTBP_Inplane_Lagrange(xL1_em, mu_em);

% Central eigenspace
eigenvec_L1_3 = eigenvec_L1(:, 3);
eigenvec_L1_4 = eigenvec_L1(:, 4);

% V real and complexe
Vr_L1 = real(eigenvec_L1_3);
Vi_L1 = abs(imag(eigenvec_L1_3));
Vr_L1_norm = Vr_L1 / norm(real(eigenvec_L1_3));
Vi_L1_norm = Vi_L1 / norm(imag(eigenvec_L1_3));

% Miscellaneous
initial_var_norm = initial_var / norm(initial_var);
initial_var_norm = initial_var_norm([1, 2, 4, 5]);

% Check
if (norm(initial_var_norm - Vr_L1_norm') < 1e-3) || (norm(initial_var_norm + Vr_L1_norm') < 1e-3)
    fprintf("Identical results for real part of the complex eigenvector.");
elseif (norm(initial_var - Vi_L1_norm) < 1e-3) || (norm(initial_var_norm + Vi_L1_norm) < 1e-3)
    fprintf("Identical results for imaginary part of the complex eigenvector.");
end

%% Problem 5
clc; close;

% Eigenvectors and eigenvalues of L4, in-plane modes
[eigenvec_L4, eigenval_L4] = Eig_CRTBP_Inplane_Lagrange(xL4_em, mu_em);
A = Jacobian_CRTBP(xL4_em, mu_em);

% Get modes
lambda1 = eigenval_L4(1,1);
V1 = eigenvec_L4(:,1);
lambda3 = eigenval_L4(3,3);
V3 = eigenvec_L4(:,3);

% Long period and short period variations
if abs(imag(lambda1))>abs(imag(lambda3))
    % Lambda 1 has the shortest period
    Vr_short = real(V1);
    Vi_short= abs(imag(V1));
    P_short = 2*pi / abs(imag(lambda1));
    Vr_long = real(V3);
    Vi_long = abs(imag(V3));
    P_long = 2*pi / abs(imag(lambda3));
else
    % Lambda 1 has the shortest period
    Vr_short = real(V3);
    Vi_short= abs(imag(V3));
    P_short = 2*pi / abs(imag(lambda3));
    Vr_long = real(V1);
    Vi_long = abs(imag(V1));
    P_long = 2*pi / abs(imag(lambda1));
end

% Scaling and enlarging
scale = 0.02;
Vr_short_scaled_4 = ScaleVariation(Vr_short, scale);
Vr_long_scaled_4 = ScaleVariation(Vr_long, scale);
Vi_short_scaled_4 = ScaleVariation(Vi_short, scale);
Vi_long_scaled_4 = ScaleVariation(Vi_long, scale);
Vr_short_scaled = [Vr_short_scaled_4(1:2)', 0, Vr_short_scaled_4(3:4)', 0];
Vr_long_scaled = [Vr_long_scaled_4(1:2)', 0, Vr_long_scaled_4(3:4)', 0];
Vi_short_scaled = [Vi_short_scaled_4(1:2)', 0, Vi_short_scaled_4(3:4)', 0];
Vi_long_scaled = [Vi_long_scaled_4(1:2)', 0, Vi_long_scaled_4(3:4)', 0];

% Short period with Complex
% Integrate variation CR3BP and obtain trajectory
options = odeset('RelTol', 2.22045*1e-14, 'AbsTol', 2.22045*1e-16);
[~, xx_short_var] = ode113(@(t,X) CRTBP_lin(t, X, A), [0 P_short],...
    Vr_short_scaled, options);
xx_short_lin = xx_short_var;
for j = 1:length(xx_short_var)
    xx_short_lin(j, 1:3) = xx_short_lin(j, 1:3) + xL4_em;
end
% Integrate CR3BP
[~, xx_short] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_short],...
    Vr_short_scaled + [xL4_em, 0, 0, 0], options);

% Plot
gca = figure(1);
plot3(xx_short_lin(:,1), xx_short_lin(:,2), xx_short_lin(:,3), 'b', 'LineWidth', 1);
hold on;
plot3(xx_short(:,1), xx_short(:,2), xx_short(:,3), 'r', 'LineWidth', 1);
plot3(xL4_em(1), xL4_em(2), xL4_em(3), 'c.', 'MarkerSize', 10);
plot3(xx_short(1,1), xx_short(1,2), xx_short(1,3), 'k*', 'Markersize', 7);
quiver3(xx_short(1,1)-0.0006, xx_short(1,2), xx_short(1,3),...
    xx_short(1,4)/2, xx_short(1,5)/2, xx_short(1,6)/2, 'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
legend('Linearized CR3BP', 'CR3BP', '$L_4$', '$\mathbf{x}_0$',...
    'interpreter','latex',...
    'Location','southwest', 'FontSize', 10);
grid on;
real axis;
view(2);

% Save plot
exportgraphics(gca, 'Prob5short.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Long period with Complex
% Integrate variation CR3BP and obtain trajectory
options = odeset('RelTol', 2.22045*1e-14, 'AbsTol', 2.22045*1e-16);
[~, xx_long_var] = ode113(@(t,X) CRTBP_lin(t, X, A), [0 P_long],...
    Vr_long_scaled, options);
xx_long_lin = xx_long_var;
for j = 1:length(xx_long_var)
    xx_long_lin(j, 1:3) = xx_long_lin(j, 1:3) + xL4_em;
end
% Integrate CR3BP
[~, xx_long] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_long],...
    Vr_long_scaled + [xL4_em, 0, 0, 0], options);

% Plot
gca = figure(2);
plot3(xx_long_lin(:,1), xx_long_lin(:,2), xx_long_lin(:,3), 'b', 'LineWidth', 1);
hold on;
plot3(xx_long(:,1), xx_long(:,2), xx_long(:,3), 'r', 'LineWidth', 1);
plot3(xL4_em(1), xL4_em(2), xL4_em(3), 'c.', 'MarkerSize', 10);
plot3(xx_long(1,1), xx_long(1,2), xx_long(1,3), 'k*', 'Markersize', 7);
quiver3(xx_long(1,1)+0.0006, xx_long(1,2), xx_long(1,3),...
    3*xx_long(1,4), 3*xx_long(1,5), 3*xx_long(1,6), 'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
legend('Linearized CR3BP', 'CR3BP', '$L_4$', '$\mathbf{x}_0$',...
    'interpreter','latex',...
    'Location','southwest', 'FontSize', 10);
grid on;
real axis;
view(2);

% Save plot
exportgraphics(gca, 'Prob5long.pdf', 'ContentType','image',...
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


% Compute Eigenvalues/Eigenvectors In-plane modes
function [eigenvectors, eigenvalues] = Eig_CRTBP_Inplane_Lagrange(X, mu)

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
df2dy = 1-(1-mu)/r1_norm^3+3*(1-mu)*y^2/r1_norm^5-mu/r2_norm^3+...
    3*mu*y^2/r2_norm^5;

% Jacobian
A_inplane = [0, 0, 1, 0;...
    0, 0, 0, 1;...
    df1dx, df1dy, 0, 2;...
    df1dy, df2dy, -2, 0];


% Eigenvalues
[eigenvectors, eigenvalues] = eig(A_inplane);

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


% Compute analytically initial guess for Periodic Orbits about L1
function [initial_var, initial_guess, xL1, P, A] =...
    init_cond_periodic_planarL1(xi0, eta0, mu)

% Compute equilibrium point, Jacobian and eigenvalues
[xL1, ~, ~, ~, ~] = EquilibriumPoints(mu);
A = Jacobian_CRTBP(xL1, mu);
[eig_L1_inplane, ~] = Eigenvalues_CRTBP_Lagrange(xL1, mu);

% Miscellaneous
lam3sq = (eig_L1_inplane(3))^2;
Uxx_eq = A(4,1);
P = 2*pi / abs(imag(eig_L1_inplane(3)));

% Compute initial variations
xidot0 = 2 * lam3sq * eta0 / (lam3sq - Uxx_eq);
etadot0 = 0.5 * (lam3sq - Uxx_eq) * xi0;
initial_var = [xi0, eta0, 0, xidot0, etadot0, 0];

% Compute non-linear initial guess
initial_guess = [xL1(1), xL1(2), xL1(3), 0, 0, 0] + initial_var;

end


% Equations of Motions CR3BP 
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

function dXdt = CRTBP_lin(~, X, A)

%Initialize
dXdt=zeros(6,1);

xi = X(1);
eta = X(2);
delta = X(3);
xidot = X(4);
etadot = X(5);
deltadot = X(6);

% Dynamics
dXdt(1:3) = [xidot; etadot; deltadot];
dXdt(4:6) = [2*etadot + A(4,1)*xi + A(4,2)*eta + A(4,3)*delta;...
    -2*xidot +  A(5,1)*xi + A(5,2)*eta + A(5,3)*delta;...
    A(6,1)*xi + A(6,2)*eta + A(6,3)*delta];

end

% Scale vector wrt position variation
function Vscaled = ScaleVariation(V, scale)
Vscaled = scale * V / norm(V(1:2));
end

% Compute Jacobi Constant in the xy-plane
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

