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

%% Exercise 1
clc; close;

%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[xL1, ~, ~, ~, ~] = EquilibriumPoints(mu_em);
x0guess_L1 = [0.836865132366261; 0; 0; 0; 0.000418613659743; 0];
P0guess_L1 = 2.691579560801770;

% Random integration with STM
[tt, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), [0 P0guess_L1],...
[x0guess_L1; reshape(eye(6),[],1)], options);

%Plot trajectory
gca = figure(1);
plot3(xx(:,1), xx(:,2),xx(:,3),'b','LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
axis([0.9999 * min(xx(:,1)), 1.0001 * max(xx(:,1)), 1.1 * min(xx(:,2)),...
    1.1 * max(xx(:,2)), -1, 1]);
view(2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xx(1,1), xx(1,2),xx(1,3),'*k','MarkerSize', 6);
quiver3(xx(1,1)-0.000008, xx(1,2), xx(1,3),...
    xx(1,4)/4, xx(1,5)/4, xx(1,6)/4,...
    'k', 'LineWidth', 1.5);
legend('Trajectory', '$L_1$', '$\mathbf{x}_0$', 'interpreter', 'latex',...
    'Location', 'northeast');
view(2);

% Export graphics
exportgraphics(gca, 'Prob1A.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot STM
gca = figure(2);
count = floor(length(tt) / 3);

% Generate state transition matrix plot for first time step
STM_iter = (reshape(xx(1,7:42),6,6))';
subplot(1, 2, 1);
imagesc(STM_iter);
colorbar;             
title(['Time Step: ', num2str(1)], 'Interpreter', 'latex');

% Generate state transition matrix plot for final time step
STM_iter = (reshape(xx(end, 7:42), 6, 6))';
subplot(1, 2, 2);
imagesc(STM_iter);
colorbar;             % Add colorbar to show scale
title(['Time Step: ', num2str(length(tt))], 'Interpreter', 'latex');

% Export graphics
exportgraphics(gca, 'Prob1B.pdf', 'ContentType','image',...
    'Resolution', 1000);


%% Exercise 2
clc; close;

% Correction scheme
tol = 1e-10;
[x0_corr, P_corr, normF_iter] = general_corrector(x0guess_L1,...
    P0guess_L1, mu_em, tol);

%Final integration
[~, xx_corr] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 P_corr],[x0_corr; reshape(eye(6),[],1)], options);

%Plot trajectory
gca = figure(1);
plot3(xx_corr(:,1), xx_corr(:,2),xx_corr(:,3),'b','LineWidth', 1.5);
hold on; 
plot3(xx(:,1), xx(:,2),xx(:,3),'r','LineWidth',1.5);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0.5 0.5 1]);
plot3(xx(1,1), xx(1,2),xx(1,3),'*k','MarkerSize', 6);
quiver3(xx(1,1)-0.000008, xx(1,2), xx(1,3),...
    xx(1,4)/4, xx(1,5)/4, xx(1,6)/4,...
    'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
axis([0.9999 * min(xx(:,1)), 1.0001 * max(xx(:,1)), 1.1 * min(xx(:,2)),...
    1.1 * max(xx(:,2)), -1, 1]);
legend('Corrected', 'Guess', '$L_1$', '$\mathbf{x}_0$',...
    'Interpreter', 'latex', 'Location','northeast');
view(2);

% Export graphics
exportgraphics(gca, 'Prob2C.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot correction performance
iter = 1:length(normF_iter);
gca = figure(2);
semilogy(iter, normF_iter, 'k*-', 'LineWidth', 1.5);
hold on;
semilogy(iter, tol * ones(length(normF_iter), 1), 'r--',...
    'LineWidth', 1.5);
xlabel('Iterations [-]', 'interpreter','latex');
ylabel('$|\mathbf{F}_k|$ [-]', 'interpreter','latex');
grid on;
real axis;
xlim([1, length(normF_iter)]);
legend('Norm constraint vector', 'Tolerance', 'Interpreter', 'latex');

% Export graphics
exportgraphics(gca, 'Prob2D.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 3
clc; close;

% Continuation scheme (Input: already corrected orbit FROM PROBLEM 2)
tol = 1e-10;  
step = 1e-4;
iters = 1000; 
x_fin = x0_corr(1) - step * iters; 
[X0_vec, P_vec] = xnatural_continuation(x0_corr, P_corr,...
    x0_corr(1), x_fin, iters, mu_em, tol, true);

% Integrate and plot 
disc_plot = 10; 
indices = round(linspace(1, iters, disc_plot));
colormap_vals = linspace(0, 1, disc_plot);
colormap_custom = [colormap_vals.', zeros(disc_plot, 1),...
    1 - colormap_vals.'];
gca1 = figure(1);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
hold on;
for i = indices
    %Final integration
    [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 P_vec(i)],[X0_vec(:,i); reshape(eye(6),[],1)], options);
    % Plot
    norm_x0 = (X0_vec(1, i) - min(X0_vec(1, :))) /...
    (max(X0_vec(1, :)) - min(X0_vec(1, :)));
    color_val = norm_x0;
    plot3(xx(:,1), xx(:,2), xx(:,3), 'Color',...
        colormap_custom(indices == i,:), 'LineWidth',1.5);
end
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
c = colorbar('eastoutside');
c.Label.String = 'Continuated $x_0$';
c.Label.Interpreter = 'latex';
colormap(c, colormap_custom);
c.TickLabelInterpreter = 'latex'; 
c.Ticks = linspace(0, 1, disc_plot);
c.TickLabels = cellstr(num2str(X0_vec(1, indices).'));
grid on;
axis equal; 
legend('$L_1$', 'Interpreter', 'latex', 'Location', 'northeast');
axis([0.99 * min(xx(:,1)), 1.01 * max(xx(:,1)), 1.1 * min(xx(:,2)),...
    1.1 * max(xx(:,2)), -1, 1]);
view(2);

% Export graphics
exportgraphics(gca1, 'Prob3A.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot: Period - Jacobi Constants
adim_to_days = 3.751902619518436 * 1e5 / 86400;
C_vec = zeros(iters, 1);
for i = 1:iters
    C_vec(i) = JacobiConstant(X0_vec(:,i), mu_em);
end
gca2 = figure(2);
plot(P_vec * adim_to_days, C_vec, 'r', 'LineWidth', 1.5);
hold on;
xlabel('$P$ [days]', 'interpreter','latex');
ylabel('$C_J$ [-]', 'interpreter','latex');
real axis;
grid on;

% Export graphics 
exportgraphics(gca2, 'Prob3B.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 4
clc; close;

% Correction scheme
tol = 1e-10;
x0guess = [0.82340, 0, -0.026755, 0, 0.13742,0]';
P0guess = 2.7477;
[x0_corr, P_corr, normF_iter] = general_corrector(x0guess,...
    P0guess, mu_em, tol);

%Final integrations
[~, xx_guess] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 P0guess],[x0guess; reshape(eye(6),[],1)], options);
[~, xx_corr] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
    [0 P_corr],[x0_corr; reshape(eye(6),[],1)], options);

%Plot trajectory
gca = figure(1);
plot3(xx_corr(:,1), xx_corr(:,2),xx_corr(:,3),'b','LineWidth', 1.5);
hold on; 
plot3(xx_guess(:,1), xx_guess(:,2),xx_guess(:,3),'r','LineWidth',1.5);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0.5 0.5 1]);
plot3(xx_guess(1,1), xx_guess(1,2),xx_guess(1,3),'*k','MarkerSize', 6);
quiver3(xx_guess(1,1)-0.008, xx_guess(1,2), xx_guess(1,3),...
    xx_guess(1,4)/4, xx_guess(1,5)/4, xx_guess(1,6)/4,...
    'k', 'LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
% axis([0.9 * min(xx_guess(:,1)), 1.1 * max(xx_guess(:,1)),...
%     1.1 * min(xx_guess(:,2)), 1.1 * max(xx_guess(:,2)), -1, 1]);
legend('Corrected', 'Guess', '$L_1$', '$\mathbf{x}_0$',...
    'Interpreter', 'latex', 'Location','northeast');

% Export graphics
exportgraphics(gca, 'Prob4A.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot correction performance
iter = 1:length(normF_iter);
gca = figure(2);
semilogy(iter, normF_iter, 'k*-', 'LineWidth', 1.5);
hold on;
semilogy(iter, tol * ones(length(normF_iter), 1), 'r--',...
    'LineWidth', 1.5);
xlabel('Iterations [-]', 'interpreter','latex');
ylabel('$|\mathbf{F}_k|$ [-]', 'interpreter','latex');
grid on;
real axis;
xlim([1, length(normF_iter)]);
legend('Norm constraint vector', 'Tolerance', 'Interpreter', 'latex');

% Export graphics
exportgraphics(gca, 'Prob4B.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Problem 5
clc; close;

% Continuation scheme (Input: already corrected orbit FROM PROBLEM 4)
tol = 1e-10;  
step = 5 * 1e-3;  
iters = 500; 
[X0_vec, P_vec] = xpseudoarc_continuation(x0_corr, P_corr,...  
      step, iters, mu_em, tol, 1, true); 

% Integrate and plot 
disc_plot = 30; 
indices = round(linspace(1, iters, disc_plot));
colormap_vals = linspace(0, 1, disc_plot);
colormap_custom = [colormap_vals.', zeros(disc_plot, 1),...
    1 - colormap_vals.'];
gca1 = figure(1);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
hold on;
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
legend('$L_1$', 'Interpreter', 'latex', 'Location', 'northeast');

% Export graphics
exportgraphics(gca1, 'Prob5A.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Plot: Period - Jacobi Constants
adim_to_days = 3.751902619518436 * 1e5 / 86400;
C_vec = zeros(iters, 1);
for i = 1:iters
    C_vec(i) = JacobiConstant(X0_vec(:,i), mu_em);
end
gca2 = figure(2);
plot(P_vec * adim_to_days, C_vec, 'r', 'LineWidth', 1.5);
hold on;
xlabel('$P$ [days]', 'interpreter','latex');
ylabel('$C_J$ [-]', 'interpreter','latex');
real axis;
grid on;

% Export graphics 
exportgraphics(gca2, 'Prob5B.pdf', 'ContentType','image',...
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

% Natural parameter continuation on x0 parameter (modified vector)
function [X0_vec, P_vec] = xnatural_continuation(X0_in, P_in, x_in,...
    x_fin, disc, mu, tol, print)
% Initialization
if nargin < 8
    print = false;
end
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
x_vec = linspace(x_in, x_fin, disc);
X0_vec = zeros(6, disc);
P_vec = zeros(disc, 1);
iter_cont = 1;
V = [X0_in; P_in];
X0_vec(:, 1) = V(1:6);  % OSS: here, change only P and not V a-priori.
P_vec(1) = V(7);
mat_mod = [eye(4), zeros(4,2); zeros(1,5), 1];

for x0 = x_vec(2:end)
    % Initialization first loop
    H = ones(7,1);
    iter_corr = 1;
    flag = 0;

    while norm(H) > tol
        %Integration state and STM
        [~, xx] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu), [0 V(end)],...
            [V(1:6); reshape(eye(6), [], 1)], options);
    
        %Reshape final STM
        STMf=(reshape(xx(end, 7:42), 6, 6))';
    
        %Jacobian of H
        xMdotf = STMexact_CRTBP(V(end), xx(end, :), mu);
        DH = [STMf([1:4, 6], :) - mat_mod, xMdotf([1:4, 6]);...
            0, 1, 0, 0, 0, 0, 0; ...
            1, 0, 0, 0, 0, 0, 0];
    
        % Update H
        H(1:5) = xx(end, [1:4, 6])- xx(1, [1:4, 6]);
        H(6) = xx(1, 2);
        H(7) = V(1) - x0;
        
        % Not modified formulation
        % DH = [STMf - eye(6), xMdotf(1:6);...
        %     1, 0, 0, 0, 0, 0, 0];
        % 
        % % Update H
        % H(1:6) = xx(end,1:6)- xx(1,1:6);
        % H(7) = V(1) - x0;

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

    % Final correction
    iter_cont = iter_cont + 1;
    X0_vec(:, iter_cont) = V(1:6);
    P_vec(iter_cont) = V(7);
end
end

% Pseudo-Arclength continuation on x0 parameter (modified vector)
function [X0_vec, P_vec] = xpseudoarc_continuation(X0_in, P_in, ...
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
initial_guess = [xL1(1), xL1(2), xL1(3), 0, 0, 0] + initial_var;

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



