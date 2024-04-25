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

% Initialization
disc = 300;
P_map = 600;
Cj = 3.175;
[xL1, xL2, xL3, xL4, xL5] = EquilibriumPoints(mu_em);
X0 = Poincare1_ICs(mu_em, Cj, disc);
options_event = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events', @event);

% Check Jacobi Constant
Cj_vec = zeros(disc, 1);
for i = 1:disc
    Cj_vec(i) = JacobiConstant(X0(:,i), mu_em);
end

% Plot trajectories
gca1 = figure(1);
ZVCxy(mu_em, Cj);
hold on; 
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 's', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
plot3(ones(100,1) * (1-mu_em), linspace(-0.12, 0.12, 100),...
    zeros(100,1), 'k--', 'LineWidth', 1.5);
plot3(X0(1,:), X0(2,:), X0(3,:), 'xk', 'MarkerSize', 8);
quiver3(X0(1, disc/4) + 0.01, X0(2, disc/4), X0(3, disc/4),...
    X0(4, disc/4) / 8, X0(5, disc/4) / 8, X0(6, disc/4) / 8,...
    'k', 'LineWidth', 1.5);
quiver3(X0(1, 3*disc/4) + 0.01, X0(2, 3*disc/4), X0(3, 3*disc/4),...
    X0(4, 3*disc/4) / 8, X0(5, 3*disc/4) / 8, X0(6, 3*disc/4) / 8,...
    'k', 'LineWidth', 1.5);
xlim([xL1(1), xL2(1)]);
ylim([-0.12, 0.12]);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
grid on;
real axis; 
legend('ZVC', '$L_1$', '$L_2$', 'Moon', '$\Sigma$', 'ICs',...
    'Interpreter', 'latex', 'Location', 'northeast');
view(2);

% Export graphics
exportgraphics(gca1, 'Prob1a1.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Create Poincare' map
mk_size = 1;
num_crossing = zeros(disc, 1);
[~, ~, te, xe, ie] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_map], ...
    X0(:, 1), options_event);
num_crossing(1) = length(ie);
[~, ~, te2, xe2, ie2] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_map], ...
    X0(:, disc/2+1), options_event);
gca2 = figure(2);
plot(xe(2:end,5), xe(2:end,2), '.b', 'MarkerSize', mk_size);
hold on;
plot(xe2(2:end,5), xe2(2:end,2), '.r', 'MarkerSize', mk_size);
for i=2:disc
    [~, ~, te, xe, ie] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_map], ...
        X0(:,i), options_event);
    if i < (disc/2 + 1)
        plot(xe(2:end, 5), xe(2:end, 2), '.b', 'MarkerSize', mk_size);
    else
        plot(xe(2:end, 5), xe(2:end, 2), '.r', 'MarkerSize', mk_size);
    end
    num_crossing(i) = length(ie);
    disp(i);
end
xlabel('$\dot{y}$ [-]', 'interpreter','latex');
ylabel('$y$ [-]', 'interpreter','latex');
real axis; 
xlim([-2.5, 2.5]);
ylim([-0.12, 0.12]);
legend('Retrograde IC', 'Prograde IC','Interpreter', 'latex',...
    'Location', 'northeast');

% Export graphics
exportgraphics(gca2, 'Prob1a2.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Analysis number of crossing
gca3 = figure(3);
plot(X0(2, :), num_crossing, '.k', 'MarkerSize', 8);
xlabel('$y_0$ [-]', 'interpreter','latex');
ylabel('Number of Crossings [-]', 'interpreter','latex');
xlim([-0.12, 0.12]);
real axis;
grid on;

% Export graphics
exportgraphics(gca3, 'Prob1a3.pdf', 'ContentType','image',...
    'Resolution', 1000);

%% Exercise 2
clc; close;

% Initial conditions
x0_1 = [0.8213849, 0, 0, 0, 0.1475143, 0];
x0_2 = [1.164855, 0, 0, 0, -0.0516671, 0];
P_1 = 2.763299;
P_2 = 3.377214;

% Correction scheme
tol = 1e-10;
[x0_1, P_1, normF_iter_1] = general_corrector(x0_1',...
    P_1, mu_em, tol);
[x0_2, P_2, normF_iter_2] = general_corrector(x0_2',...
    P_2, mu_em, tol);

% Check Jacobi Constant
C_1 = JacobiConstant(x0_1, mu_em);
C_2 = JacobiConstant(x0_2, mu_em);


% Plot orbits
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[tt1, xx1] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_1], x0_1, options);
[tt2, xx2] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_2], x0_2, options);

%Plot trajectory
gca = figure(1);
ZVCxy(mu_em, C_2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(xx1(:,1), xx1(:,2), xx1(:,3),'b','LineWidth', 1.5);
plot3(xx2(:,1), xx2(:,2), xx2(:,3),'r','LineWidth', 1.5);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
xlim([0.8, 1.2]);
ylim([-0.1, 0.1]);
legend('ZVC', '$L_1$', '$L_2$', 'interpreter', 'latex',...
    'Location', 'northeast');
view(2);

% Export graphics
exportgraphics(gca, 'Prob2a.pdf', 'ContentType','image',...
    'Resolution', 1000);

% Generate initial conditions manifolds
options_Moon_for = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events', @event_Moon_for);
options_Moon_back = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events', @event_Moon_back);
l_star_em = 384400;
d = 0.1 / l_star_em; 
disc = 300;
T_man = 10;
[~, X0_man_unstab_pos_1] = ic_manifolds_pos(x0_1, P_1, d, disc,...
    mu_em);
[X0_man_stab_neg_2, ~] = ic_manifolds_neg(x0_2, P_2, d, disc,...
    mu_em);

% Initialize to record crossings
num_crossing = zeros(2, disc-1);

% Plot Poincarè
gca2 = figure(2);
for i=2:disc
    % Integration with STM - positive manifolds
    [~, xx_unstab_pos_1, te1, xe1, ie1] = ode113(@(t, X) CRTBP(t, X, mu_em),...
        [0 T_man], X0_man_unstab_pos_1(:, i), options_Moon_for);
    [~, xx_stab_pos_2, te2, xe2, ie2] = ode113(@(t, X) CRTBP(t, X, mu_em),...
        [0 -T_man], X0_man_stab_neg_2(:, i), options_Moon_back);

   % Retain only the first two elements 
   te1 = te1(1:min(2, length(ie1)));
   xe1 = xe1(1:min(2, length(ie1)),:);
   ie1 = ie1(1:min(2, length(ie1)));
   te2 = te2(1:min(2, length(ie2)));
   xe2 = xe2(1:min(2, length(ie2)),:);
   ie2 = ie2(1:min(2, length(ie2)));

   % Record number of crossings
   num_crossing(1, i-1) = length(ie1);
   num_crossing(2, i-1) = length(ie2);

    % Plot
    plot(xe1(:,5), xe1(:, 2), '.r', 'MarkerSize', 6);
    hold on;
    plot(xe2(:,5), xe2(:,2), '.b', 'MarkerSize', 6);
end
xlabel('$\dot{y}$ [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
xlim([-0.6, 3.45]);
ylim([-0.085, 0]);
legend('$L_1$, Unstable Manifold', '$L_2$, Stable Manifold', 'interpreter', 'latex',...
    'Location', 'southeast');

% Export graphics
exportgraphics(gca2, 'Prob2b.pdf', 'ContentType','image',...
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
    % Data
    GM_earth = 398600.435436; % km^3/s^2
    GM_moon = 4902.800066; 
    mu_em = GM_moon / (GM_earth + GM_moon);
    
    % Defintion value and direction
    value = X(1) - (1-mu_em);   % x = 1-mu
    direction = 1; % Event function increasing
    isterminal = 0;
    
    % Definition isterminal
    % persistent count;
    % if isempty(count)
    %     count = 300;
    % end
    % if count > 1
    %     count = count - 1;
    %     isterminal = 0;  
    % else
    %     isterminal = 1; % Stop 
    % end
    % disp(count);
end

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

% Compute Psuedo-Potential
function U_star = PseudoPotential(X, mu)
    % Initialization
    x = X(1);
    y = X(2);
    z = X(3);
    
    % Jacobi Constant Computation
    U_star = 0.5 * (x^2+y^2) + (1-mu)/sqrt((x+mu)^2+y^2+z^2) + ...
        mu/sqrt((x-1+mu)^2+y^2+z^2);
end

% Poincare map initial conditions
function X0 = Poincare1_ICs(mu, C, disc)
    % Initialization
    X0 = zeros(6, disc);
    f_yzvc = @(y) ((1-mu)^2+y^2) + 2*(1-mu)/sqrt(1+y^2) + ...
        2*mu/sqrt(y^2) - C;
    y_zvc = fsolve(f_yzvc, 0.1);
    l_star = 384400;
    r_moon = 1737.4;
    
    % Select position components
    X0(1, :) = 1-mu; % ICs on the Sigma-plane
    X0(2, 1:end/2) = linspace(2 * r_moon / l_star, 0.95 * y_zvc, disc/2);
    X0(2, (end/2+1):end) = linspace(- 2 * r_moon / l_star, - 0.95 * y_zvc, disc/2);
    X0(3, :) = 0; % ICs on the xy-plane
    
    % Select velocity components
    % v_square = zeros(disc, 1);
    % for i = 1:disc
    %     x = X0(1, i);
    %     y = X0(2, i);
    %     z = X0(3, i);
    %     v_square(i) = (x^2+y^2) + 2*(1-mu)/sqrt((x+mu)^2+y^2+z^2) +...
    %         2*mu/sqrt((x-1+mu)^2+y^2+z^2) - C;
    % end
    for i = 1:disc
        X0(4, i) = sqrt(2 * PseudoPotential(X0(:, i), mu) - C);
    end
    X0(5, :) = 0;
    X0(6, :) = 0; 

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

% Events Moon
function [value, isterminal, direction] = event_Moon_for(~, X, ~)
% Data
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 

% Computations mass ratios
mu_em = GM_moon / (GM_earth + GM_moon);

% Event
value = X(1) - (1-mu_em);   % x = 1-mu
isterminal = 0;
direction = 1; % all directions
end

function [value, isterminal, direction] = event_Moon_back(~, X, ~)
% Data
GM_earth = 398600.435436; % km^3/s^2
GM_moon = 4902.800066; 

% Computations mass ratios
mu_em = GM_moon / (GM_earth + GM_moon);

% Event
value = X(1) - (1-mu_em);   % x = 1-mu
isterminal = 0;
direction = -1; % all directions
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

