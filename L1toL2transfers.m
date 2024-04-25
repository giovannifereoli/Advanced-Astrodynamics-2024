%% Advanced Astrodynamics (ASEN6060) - Project, First Half
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
l_star = 384400;  % km
t_star = 3.751902619518436 * 1e5; % km/s

% Computations mass ratios
mu_se = GM_earth / (GM_sun + GM_earth); % Mass Ratio
mu_em = GM_moon / (GM_earth + GM_moon);

%% Plot Lyapunov Orbits
clc; close;

% Initial conditions
[xL1, xL2, xL3, xL4, xL5] = EquilibriumPoints(mu_em);
x0_1 = [0.821950426219030, 0, 0, 0, 0.141479662833491, 0];
x0_2 = [1.175773196736922, 0, 0,  0, -0.119977116007445, 0];
P_1 = 2.757108054159905;
P_2 = 3.396688765837098;

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
[~, xx1] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_1], x0_1, options);
[~, xx2] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 P_2], x0_2, options);

%Plot trajectory
gca = figure(1);
ZVCxy(mu_em, C_2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.8 0.2]);
plot3(xx1(:,1), xx1(:,2), xx1(:,3),'k','LineWidth', 2.5);
plot3(xx2(:,1), xx2(:,2), xx2(:,3),'k','LineWidth', 2.5);
quiver3(xx1(1,1) - 0.01, xx1(1,2), xx1(1,3),...
    xx1(1,4) / 3, xx1(1,5) / 3, xx1(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(xx2(1,1) + 0.01, xx2(1,2), xx2(1,3),...
    xx2(1,4) / 3, xx2(1,5) / 3, xx2(1,6) / 3,...
    'k', 'LineWidth', 2);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
xlim([0.8, 1.2]);
ylim([-0.1, 0.1]);
legend('ZVC ($C_{J,2})$', '$L_1$', '$L_2$', 'Moon', 'interpreter', 'latex',...
    'Location', 'eastoutside', 'FontSize', 8);
view(2);

% Export graphics
% exportgraphics(gca, 'Project1_orbits.pdf', 'ContentType','image',...
%     'Resolution', 1000);

%% Plot Manifolds
clc; close;

%Initialization
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
l_star_em = 384400;
d = 0.1 / l_star_em; 
disc = 600;
T_man_1 = 5.5;
T_man_2 = 7;

% Generate initial conditions manifolds
[~, X0_man_unstab_pos_1] = ic_manifolds_pos(x0_1, P_1, d, disc,...
    mu_em);
[X0_man_stab_neg_2, ~] = ic_manifolds_neg(x0_2, P_2, d, disc,...
    mu_em);

% Integration with STM
[~, xx_unstab_pos_1] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), ...
    [0 T_man_1], [X0_man_unstab_pos_1(:, 1); reshape(eye(6),[],1)], options);
[~, xx_stab_neg_2] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em), ...
    [0 -T_man_2], [X0_man_stab_neg_2(:, 1); reshape(eye(6),[],1)], options);

%Plot trajectory
gca2 = figure(2);
plot3(xx_stab_neg_2(:,1), xx_stab_neg_2(:,2),xx_stab_neg_2(:,3),'b',...
    'LineWidth', 0.5);
hold on;
plot3(xx_unstab_pos_1(:,1), xx_unstab_pos_1(:,2),xx_unstab_pos_1(:,3),'r',...
    'LineWidth', 0.5);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.2 0.5 1]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.3 1]);
for i=2:disc
    % Integration with STM - positive manifolds
    [~, xx_unstab_pos_1] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 T_man_1], [X0_man_unstab_pos_1(:, i); reshape(eye(6),[],1)], options);
    [~, xx_stab_neg_2] = ode113(@(t,xM) STMexact_CRTBP(t, xM, mu_em),...
        [0 -T_man_2], [X0_man_stab_neg_2(:, i); reshape(eye(6),[],1)], options);

    % Plot
    plot3(xx_stab_neg_2(:,1), xx_stab_neg_2(:,2),xx_stab_neg_2(:,3),'b',...
        'LineWidth', 0.5);
    plot3(xx_unstab_pos_1(:,1), xx_unstab_pos_1(:,2),xx_unstab_pos_1(:,3),'r',...
        'LineWidth', 0.5);
    disp(i);
end
plot3(xx1(:,1), xx1(:,2), xx1(:,3),'k','LineWidth', 2.5);
plot3(xx2(:,1), xx2(:,2), xx2(:,3),'k','LineWidth', 2.5);
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.3 1]);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
view(2);
legend('Stable Manifold', 'Unstable Manifold',...
    '$L_1$', '$L_2$', 'Moon', 'interpreter', 'latex', 'Location',...
    'eastoutside', 'FontSize', 8);

% Export graphics
% exportgraphics(gca2, 'Project1_manifolds.pdf', 'ContentType','image',...
%     'Resolution', 1000);

%% Plot Poincarè map 
clc; close;

% Generate initial conditions manifolds
options_Moon_for = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events', @event_Moon_for);
options_Moon_back = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15,...
    'Events', @event_Moon_back);
T_man = 10;

% Initialize to record crossings
num_crossing = zeros(2, disc-1);

% Initialize cell arrays to store variables for each i
xe1_cell = cell(disc-1, 1);
xe2_cell = cell(disc-1, 1);
te1_cell = cell(disc-1, 1);
te2_cell = cell(disc-1, 1);
ie1_cell = cell(disc-1, 1);
ie2_cell = cell(disc-1, 1);

% Plot Poincarè
gca3 = figure(3);
for i=2:disc
   % Integration with STM - positive manifolds
   [~, xx_unstab_pos_1, te1, xe1, ie1] = ode113(@(t, X) CRTBP(t, X, mu_em),...
       [0 T_man], X0_man_unstab_pos_1(:, i), options_Moon_for);
   [~, xx_stab_neg_2, te2, xe2, ie2] = ode113(@(t, X) CRTBP(t, X, mu_em),...
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

   % Save variables for each i
   xe1_cell{i-1} = xe1;
   xe2_cell{i-1} = xe2;
   te1_cell{i-1} = te1;
   te2_cell{i-1} = te2;
   ie1_cell{i-1} = ie1;
   ie2_cell{i-1} = ie2;

   % Plot
   plot(xe1(:,5), xe1(:, 2), '.r', 'MarkerSize', 6);
   hold on;
   plot(xe2(:,5), xe2(:,2), '.b', 'MarkerSize', 6);
   disp(i);
end
xlabel('$\dot{y}$ [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
xlim([-1.2, 3.45]);
ylim([-0.095, 0]);
legend('$L_1$, Unstable Manifold', '$L_2$, Stable Manifold', 'interpreter', 'latex',...
    'Location', 'southeast');

% Export graphics
% exportgraphics(gca3, 'Project1_Poincare.pdf', 'ContentType','image',...
%     'Resolution', 1000);

%% Find Crossings Poincarè

% Initialize
intersection_Poincare = cell(disc-1, disc-1, 4);
intersection_ICs = cell(disc-1, disc-1, 2);
threshold = 1e-3; 
counter = 0;

% Fill crossings matrix
for i = 1:disc-1
    for j = 1:disc-1
        for k = 1:size(xe1_cell{i}, 1)
            for w = 1:size(xe2_cell{j}, 1)
                if norm(xe1_cell{i}(k, [2,5]) - xe2_cell{j}(w, [2,5])) < threshold
                    intersection_Poincare{i, j, 1} = xe1_cell{i}(k, :); 
                    intersection_Poincare{i, j, 2} = xe2_cell{j}(w, :);
                    intersection_Poincare{i, j, 3} = te1_cell{i}(k); 
                    intersection_Poincare{i, j, 4} = te2_cell{j}(w);
                    counter = counter +1;
                    disp(counter);
                end
            end
        end
    end
end

% Print not-zero elements
for i = 1:size(intersection_Poincare, 1)
    for j = 1:size(intersection_Poincare, 2)
            % Check if the current element is not empty
            if ~isempty(intersection_Poincare{i, j, 1})
                % Print the indices of non-empty elements
                disp(['Non-empty element found at index (', num2str(i), ',' ...
                    '', num2str(j), ')']);
                intersection_ICs{i, j, 1} = X0_man_unstab_pos_1(:, i)';
                intersection_ICs{i, j, 2} = X0_man_stab_neg_2(:, j)';
                intersection_ICs{i, j, 3} = P_1 * (i-1) / disc;
                intersection_ICs{i, j, 4} = P_2 * (j-1) / disc;
            end
    end
end

%% First Transfer

% Initialize Poincarè Intersection
row = 383; col = 379;
% row = 523; col = 176;

% Initial guesses
x0_L1 = x0_1; 
P_L1 = intersection_ICs{row, col, 3};  
x0_man_L1 = intersection_ICs{row, col, 1};
P_man_L1 = intersection_Poincare{row, col, 3};
x0_L2 = x0_2;
P_L2 = intersection_ICs{row, col, 4};
x0_man_L2 = intersection_ICs{row, col, 2};
P_man_L2 = intersection_Poincare{row, col, 4};
disc = 10000;

% Plot orbits
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[~, xx1] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, P_L1, disc), x0_L1, options);
[~, xx2] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, P_man_L1, disc), x0_man_L1, options);
[~, xx3] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, P_man_L2, disc), x0_man_L2, options); 
[~, xx4] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, (P_L2-P_2), disc), x0_L2, options);

%Plot trajectory
gca4 = figure(4);
ZVCxy(mu_em, C_2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.3 1]);
plot3(xx2(:,1), xx2(:,2), xx2(:,3),'r','LineWidth', 1.5);
plot3(xx3(:,1), xx3(:,2), xx3(:,3),'b','LineWidth', 1.5);
plot3(xx4(:,1), xx4(:,2), xx4(:,3),'k','LineWidth', 2.5);
plot3(xx1(:,1), xx1(:,2), xx1(:,3),'k','LineWidth', 2.5);
plot3(xx1(1,1), xx1(1,2), xx1(1,3), '.k', 'MarkerSize', 15);
plot3(xx1(end,1), xx1(end,2), xx1(end,3), '.k', 'MarkerSize', 15);
plot3(xx2(1,1), xx2(1,2), xx2(1,3), '.k', 'MarkerSize', 15);
plot3(xx2(end,1), xx2(end,2), xx2(end,3), '.k', 'MarkerSize', 15);
plot3(xx3(1,1), xx3(1,2), xx3(1,3), '.k', 'MarkerSize', 15);
plot3(xx3(end,1), xx3(end,2), xx3(end,3), '.k', 'MarkerSize', 15);
plot3(xx4(1,1), xx4(1,2), xx4(1,3), '.k', 'MarkerSize', 15);
plot3(xx4(end,1), xx4(end,2), xx4(end,3), '.k', 'MarkerSize', 15);
quiver3(xx1(1,1) - 0.01, xx1(1,2), xx1(1,3),...
    xx1(1,4) / 3, xx1(1,5) / 3, xx1(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(xx4(1,1) + 0.01, xx4(1,2), xx4(1,3),...
    xx4(1,4) / 3, xx4(1,5) / 3, xx4(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(0.965, -0.05, 0,...
    0.05, 0,  0,...
    'k', 'LineWidth', 2);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
xlim([0.8, 1.2]);
ylim([-0.12, 0.12]);
legend('ZVC ($C_{J,2})$', '$L_1$', '$L_2$', 'Moon',...
    '$L_1$, Unstable Manifold', '$L_2$, Stable Manifold',...
    'interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 8);
view(2);

% Reverse few integrations due stable manifold for NEXT section
x0_man_L2 = xx3(end, :);
P_man_L2 = - P_man_L2;
x0_L2 = xx4(end, :); % extra: P_L2 = (P_2-P_L2);
xxf_des = xx4; 
xx0_des = xx1;

% Export graphics
% exportgraphics(gca4, 'Project1_FirstTransferGuess.pdf', 'ContentType','image',...
%     'Resolution', 1000);

%% First Transfer Correction

% Initialization and initial guess correction
% OSS: the second manifold has been 'reversed'!
V_guess = [x0_man_L1, P_man_L1/2, xx2(end/2,:), P_man_L1/2,...
    x0_man_L2, P_man_L2/2,  xx3(end/2, :), P_man_L2/2];  
X0_des = xx1(end, :)';
Xf_des = x0_L2';
tol = 1e-10; 

% Correction
[V_corr, normF_iter] = transfer_ms(V_guess, X0_des, Xf_des, mu_em, tol);

% Integration orbits
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[~, xx1_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(7)], V_corr(1:6), options);
[~, xx2_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(14)], V_corr(8:13), options);
[~, xx3_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(21)], V_corr(15:20), options);
[~, xx4_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(28)], V_corr(22:27), options);

% Plot orbits
gca5 = figure(5);
ZVCxy(mu_em, C_2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.3 1]);
plot3(xx1_corr(1,1), xx1_corr(1,2), xx1_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx1_corr(:,1), xx1_corr(:,2), xx1_corr(:,3),'g','LineWidth', 2);
plot3(xx1_corr(end,1), xx1_corr(end,2), xx1_corr(end,3), '.k', 'MarkerSize', 15);
plot3(xx2_corr(1,1), xx2_corr(1,2), xx2_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx2_corr(:,1), xx2_corr(:,2), xx2_corr(:,3),'y','LineWidth', 1.5);
plot3(xx2_corr(end,1), xx2_corr(end,2), xx2_corr(end,3), '.k', 'MarkerSize', 15);
plot3(xx3_corr(1,1), xx3_corr(1,2), xx3_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx3_corr(:,1), xx3_corr(:,2), xx3_corr(:,3),'m','LineWidth', 1.5);
plot3(xx3_corr(end,1), xx3_corr(end,2), xx3_corr(end,3), '.k', 'MarkerSize', 15);
plot3(xx4_corr(1,1), xx4_corr(1,2), xx4_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx4_corr(:,1), xx4_corr(:,2), xx4_corr(:,3),'c','LineWidth', 2);
plot3(xx4_corr(end,1), xx4_corr(end,2), xx4_corr(end,3), '.k', 'MarkerSize', 15);
quiver3(xx1(1,1) - 0.01, xx1(1,2), xx1(1,3),...
    xx1(1,4) / 3, xx1(1,5) / 3, xx1(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(xx4(1,1) + 0.01, xx4(1,2), xx4(1,3),...
    xx4(1,4) / 3, xx4(1,5) / 3, xx4(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(0.965, -0.05, 0,...
    0.05, 0,  0,...
    'k', 'LineWidth', 2);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
xlim([0.8, 1.2]);
ylim([-0.12, 0.12]);
legend('ZVC ($C_{J,2})$', '$L_1$', '$L_2$', 'Moon', '$\Delta V$ Node', '$1^{st}$ Interval',...
    '', '', '$2^{nd}$ Interval', '', '', '$3^{rd}$ Interval', '', '', '$4^{th}$ Interval',...
    'interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 8);
view(2);

% Export graphics
% exportgraphics(gca5, 'Project1_FirstTransferCorrected.pdf', 'ContentType','image',...
%     'Resolution', 1000);

% Plot correction performance
iter = 1:length(normF_iter);
gca6 = figure(6);
semilogy(iter, normF_iter, 'k.', 'LineWidth', 0.2);
hold on;
semilogy(iter, tol * ones(length(normF_iter), 1), 'r--',...
    'LineWidth', 0.1);
xlabel('Iterations [-]', 'interpreter','latex');
ylabel('$|\mathbf{F}_k|$ [-]', 'interpreter','latex');
grid on;
real axis;
xlim([1, length(normF_iter)]);
ylim([0.1 * normF_iter(end), normF_iter(1)]);
legend('Norm constraint vector', 'Tolerance', 'Interpreter', 'latex');

% Export graphics
% exportgraphics(gca5, 'Project1_FirstTransferCorrection.pdf', 'ContentType','image',...
%     'Resolution', 1000);

% After plotting you have the flows, compute TCMs
TCM1 = xx1_corr(1, 4:6) - xx0_des(end, 4:6);
TCM2 = xx2_corr(1, 4:6) - xx1_corr(end, 4:6);
TCM3 = xx3_corr(1, 4:6) - xx2_corr(end, 4:6);
TCM4 = xx3_corr(end, 4:6) - xx4_corr(1, 4:6);
TCM5 = xx4_corr(end, 4:6) - xxf_des(1, 4:6);
TCMtot_norm = norm(TCM1) + norm(TCM2) + norm(TCM3) + norm(TCM4) + norm(TCM5);
fprintf('The trajectory TCMs need a ΔV of %s [km/s].\n',...
    num2str(TCMtot_norm * t_star / l_star));

%% Second Transfer

% Initialization
% row = 383; col = 379;
row = 523; col = 176;

% Initial guesses
x0_L1 = x0_1; 
P_L1 = intersection_ICs{row, col, 3};  
x0_man_L1 = intersection_ICs{row, col, 1};
P_man_L1 = intersection_Poincare{row, col, 3};
x0_L2 = x0_2;
P_L2 = intersection_ICs{row, col, 4};
x0_man_L2 = intersection_ICs{row, col, 2};
P_man_L2 = intersection_Poincare{row, col, 4};
disc = 10000;

% Plot orbits
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[~, xx1] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, P_L1, disc), x0_L1, options);
[~, xx2] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, P_man_L1, disc), x0_man_L1, options);
[~, xx3] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, P_man_L2, disc), x0_man_L2, options); 
[~, xx4] = ode113(@(t,X) CRTBP(t, X, mu_em), linspace(0, (P_L2-P_2), disc), x0_L2, options);

%Plot trajectory
gca7 = figure(7);
ZVCxy(mu_em, C_2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.3 1]);
plot3(xx2(:,1), xx2(:,2), xx2(:,3),'r','LineWidth', 1.5);
plot3(xx3(:,1), xx3(:,2), xx3(:,3),'b','LineWidth', 1.5);
plot3(xx4(:,1), xx4(:,2), xx4(:,3),'k','LineWidth', 2.5);
plot3(xx1(:,1), xx1(:,2), xx1(:,3),'k','LineWidth', 2.5);
plot3(xx1(1,1), xx1(1,2), xx1(1,3), '.k', 'MarkerSize', 15);
plot3(xx1(end,1), xx1(end,2), xx1(end,3), '.k', 'MarkerSize', 15);
plot3(xx2(1,1), xx2(1,2), xx2(1,3), '.k', 'MarkerSize', 15);
plot3(xx2(end,1), xx2(end,2), xx2(end,3), '.k', 'MarkerSize', 15);
plot3(xx3(1,1), xx3(1,2), xx3(1,3), '.k', 'MarkerSize', 15);
plot3(xx3(end,1), xx3(end,2), xx3(end,3), '.k', 'MarkerSize', 15);
plot3(xx4(1,1), xx4(1,2), xx4(1,3), '.k', 'MarkerSize', 15);
plot3(xx4(end,1), xx4(end,2), xx4(end,3), '.k', 'MarkerSize', 15);
quiver3(xx1(1,1) - 0.01, xx1(1,2), xx1(1,3),...
    xx1(1,4) / 3, xx1(1,5) / 3, xx1(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(xx4(1,1) + 0.01, xx4(1,2), xx4(1,3),...
    xx4(1,4) / 3, xx4(1,5) / 3, xx4(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(0.965, -0.096, 0,...
    0.05, 0,  0,...
    'k', 'LineWidth', 2);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
xlim([0.8, 1.2]);
ylim([-0.12, 0.12]);
legend('ZVC ($C_{J,2})$', '$L_1$', '$L_2$', 'Moon',...
    '$L_1$, Unstable Manifold', '$L_2$, Stable Manifold',...
    'interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 8);
view(2);

% Reverse few integrations due stable manifold for NEXT section
x0_man_L2 = xx3(end, :);
P_man_L2 = - P_man_L2;
x0_L2 = xx4(end, :);
P_L2 = (P_2-P_L2);
xxf_des = xx4; 
xx0_des = xx1;

% Export graphics
% exportgraphics(gca7, 'Project1_SecondTransferGuess.pdf', 'ContentType','image',...
%     'Resolution', 1000);

%% Second Transfer Correction

% Initialization and initial guess correction
% OSS: the second manifold has been 'reversed'!
V_guess = [x0_man_L1, P_man_L1/2, xx2(end/2,:), P_man_L1/2,...
    x0_man_L2, P_man_L2/2,  xx3(end/2, :), P_man_L2/2];
X0_des = xx1(end, :)';
Xf_des = x0_L2';
tol = 1e-10; 

% Correction
[V_corr, normF_iter] = transfer_ms(V_guess, X0_des, Xf_des, mu_em, tol);

% Integration orbits
options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-15);
[~, xx1_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(7)], V_corr(1:6), options);
[~, xx2_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(14)], V_corr(8:13), options);
[~, xx3_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(21)], V_corr(15:20), options);
[~, xx4_corr] = ode113(@(t,X) CRTBP(t, X, mu_em), [0 V_corr(28)], V_corr(22:27), options);

% Plot orbits
gca8 = figure(8);
ZVCxy(mu_em, C_2);
hold on;
plot3(xL1(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.5 0.5 1]);
plot3(xL2(1), 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.8 0.2]);
plot3(1-mu_em, 0, 0, 'd', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.3 1]);
plot3(xx1_corr(1,1), xx1_corr(1,2), xx1_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx1_corr(:,1), xx1_corr(:,2), xx1_corr(:,3),'g','LineWidth', 2);
plot3(xx1_corr(end,1), xx1_corr(end,2), xx1_corr(end,3), '.k', 'MarkerSize', 15);
plot3(xx2_corr(1,1), xx2_corr(1,2), xx2_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx2_corr(:,1), xx2_corr(:,2), xx2_corr(:,3),'y','LineWidth', 1.5);
plot3(xx2_corr(end,1), xx2_corr(end,2), xx2_corr(end,3), '.k', 'MarkerSize', 15);
plot3(xx3_corr(1,1), xx3_corr(1,2), xx3_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx3_corr(:,1), xx3_corr(:,2), xx3_corr(:,3),'m','LineWidth', 1.5);
plot3(xx3_corr(end,1), xx3_corr(end,2), xx3_corr(end,3), '.k', 'MarkerSize', 15);
plot3(xx4_corr(1,1), xx4_corr(1,2), xx4_corr(1,3), '.k', 'MarkerSize', 15);
plot3(xx4_corr(:,1), xx4_corr(:,2), xx4_corr(:,3),'c','LineWidth', 2);
plot3(xx4_corr(end,1), xx4_corr(end,2), xx4_corr(end,3), '.k', 'MarkerSize', 15);
quiver3(xx1(1,1) - 0.01, xx1(1,2), xx1(1,3),...
    xx1(1,4) / 3, xx1(1,5) / 3, xx1(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(xx4(1,1) + 0.01, xx4(1,2), xx4(1,3),...
    xx4(1,4) / 3, xx4(1,5) / 3, xx4(1,6) / 3,...
    'k', 'LineWidth', 2);
quiver3(0.965, -0.096, 0,...
    0.05, 0,  0,...
    'k', 'LineWidth', 2);
xlabel('x [-]', 'interpreter','latex');
ylabel('y [-]', 'interpreter','latex');
zlabel('z [-]', 'interpreter','latex');
grid on;
axis equal; 
xlim([0.8, 1.2]);
ylim([-0.12, 0.12]);
legend('ZVC ($C_{J,2})$', '$L_1$', '$L_2$', 'Moon', '$\Delta V$ Node', '$1^{st}$ Interval',...
    '', '', '$2^{nd}$ Interval', '', '', '$3^{rd}$ Interval', '', '', '$4^{th}$ Interval',...
    'interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 8);
view(2);

% Export graphics
% exportgraphics(gca8, 'Project1_SecondTransferCorrected.pdf', 'ContentType','image',...
%     'Resolution', 1000);

% Plot correction performance
iter = 1:length(normF_iter);
gca9 = figure(9);
semilogy(iter, normF_iter, 'k.', 'LineWidth', 0.2);
hold on;
semilogy(iter, tol * ones(length(normF_iter), 1), 'r--',...
    'LineWidth', 0.1);
xlabel('Iterations [-]', 'interpreter','latex');
ylabel('$|\mathbf{F}_k|$ [-]', 'interpreter','latex');
grid on;
real axis;
xlim([1, length(normF_iter)]);
ylim([0.1 * normF_iter(end), normF_iter(1)]);
legend('Norm constraint vector', 'Tolerance', 'Interpreter', 'latex');

% Export graphics
% exportgraphics(gca9, 'Project1_SecondTransferCorrection.pdf', 'ContentType','image',...
%     'Resolution', 1000);

% After plotting you have the flows, compute TCMs
TCM1 = xx1_corr(1, 4:6) - xx0_des(end, 4:6);
TCM2 = xx2_corr(1, 4:6) - xx1_corr(end, 4:6);
TCM3 = xx3_corr(1, 4:6) - xx2_corr(end, 4:6);
TCM4 = xx3_corr(end, 4:6) - xx4_corr(1, 4:6);
TCM5 = xx4_corr(end, 4:6) - xxf_des(1, 4:6);
TCMtot_norm = norm(TCM1) + norm(TCM2) + norm(TCM3) + norm(TCM4) + norm(TCM5);
fprintf('The trajectory TCMs need a ΔV of %s [km/s].\n',...
    num2str(TCMtot_norm * t_star / l_star));

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

% Multiple-Shooting Correction algorithm 
function [V_corr, normF_iter] = transfer_ms(V_guess, X1, X4, mu, tol)

    % Initialization
    V = V_guess;
    [F, ~] = nonlcons_ms(V, mu, X1, X4);
    normF_iter = zeros(100, 1);
    iter = 1;
    
    while norm(F) > tol
        % Jacobian of F, update F
        [F, DF] = nonlcons_ms(V, mu, X1, X4);
        normF_iter(iter) = norm(F);
        disp(normF_iter(iter))
    
        % Check if norm of F exceeds tolerance
        if norm(F) > 1e3
            fprintf('Error: Norm of F diverging. \n');
            break;  % Exit the loop if norm(F) exceeds tolerance
        end
    
        % Correction
        dV = - DF' * ((DF*DF') \ F);
        V = V + 1e-1 .* dV';
        iter = iter + 1;
        disp(iter)
    end
    
    % Final correction
    V_corr = V;
    normF_iter = normF_iter(1:iter-1);
    
    % Final print
    disp('Correction terminated successfully.');
    fprintf('         Current function value: %.14f\n', normF_iter(end));
    fprintf('         Iterations: %d\n', iter-1);
end

% Constraint vector and gradient for Multiple-Shooting
function [F, DF] = nonlcons_ms(V, mu, X0_des, Xf_des)
    %Upacking state
    X1 = V(1:6)';
    T1 = V(7);
    X2 = V(8:13)';
    T2 = V(14);
    X3 = V(15:20)';
    T3 = V(21);
    X4 = V(22:27)';
    T4 = V(28);

    %Integrations
    options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-16);

    %Flow 1
    [~, xxM1] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T1],...
            [X1; reshape(eye(6), [], 1)],...
            options);
    Xf1 = (xxM1(end, 1:6))';
    rhs_f1 = CRTBP(T1, Xf1, mu); 
    STMf1 = (reshape(xxM1(end, 7:42), 6, 6))';  
    
    %Flow 2
    [~, xxM2] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T2],...
            [X2; reshape(eye(6), [], 1)],...
            options);
    Xf2 = (xxM2(end, 1:6))';
    rhs_f2 = CRTBP(T2, Xf2, mu); 
    STMf2 = (reshape(xxM2(end, 7:42), 6, 6))';   %From equations to STM
      
    %Flow 3
    [~, xxM3] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T3],...
            [X3; reshape(eye(6), [], 1)],...
            options);
    Xf3 = (xxM3(end, 1:6))';
    rhs_f3 = CRTBP(T3, Xf3, mu); 
    STMf3 = (reshape(xxM3(end, 7:42), 6, 6))';   %From equations to STM

    %Flow 4
    [~, xxM4] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T4],...
            [X4; reshape(eye(6), [], 1)],...
            options);
    Xf4 = (xxM4(end, 1:6))';
    rhs_f4 = CRTBP(T4, Xf4, mu); 
    STMf4 = (reshape(xxM4(end, 7:42), 6, 6))';   %From equations to STM
   
    % Constraint Vector (i.e., positions on nodes)
    F = [X1(1:3) - X0_des(1:3); ...
        Xf1(1:3) - X2(1:3); ... 
        Xf2(1:3) - X3(1:3); ... 
        Xf3(1:3) - X4(1:3); ... 
        Xf4(1:3) - Xf_des(1:3)];

    %Composition of DF
    DFrow1 = [eye(3), zeros(3,3), zeros(3,1), zeros(3,6), zeros(3,1),...
        zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1)];
    DFrow2 = [STMf1(1:3, :), rhs_f1(1:3), -eye(3), zeros(3,3), zeros(3,1),...
        zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1)];
    DFrow3 = [zeros(3,6), zeros(3,1), STMf2(1:3, :), rhs_f2(1:3),...
        -eye(3), zeros(3,3), zeros(3,1), zeros(3,6), zeros(3,1)];
    DFrow4 = [zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1),...
        STMf3(1:3, :), rhs_f3(1:3), -eye(3), zeros(3,3), zeros(3,1)];
    DFrow5 = [zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1),...
        zeros(3,6), zeros(3,1), STMf4(1:3, :), rhs_f4(1:3)];

    DF = [DFrow1; DFrow2; DFrow3; DFrow4; DFrow5];

end

% EXTRA CODE to make checks (Solving NLP)
%
% options=optimoptions('fsolve', 'Algorithm', 'Levenberg-Marquardt', ...
%     'Display','iter', 'StepTolerance', 0.2 * 1e-6, ...
%     'MaxIterations', 1000, 'SpecifyObjectiveGradient', true);
% [valid,err] = checkGradients(@(V) obj_ss(V, mu_em, X0_des, Xf_des),...
%     V_guess, Display="on")
% V_corr = fsolve(@(V) obj_ss(V, mu_em, X0_des, Xf_des), V_guess, options);
% 
% function [F, DF] = obj_ss(V, mu, X0_des, Xf_des)
%    %Upacking state
%     X1 = V(1:6)';
%     T1 = V(7);
%     X2 = V(8:13)';
%     T2 = V(14);
%     X3 = V(15:20)';
%     T3 = V(21);
%     X4 = V(22:27)';
%     T4 = V(28);
% 
%     %Integrations
%     options = odeset('RelTol', 2.24*1e-14, 'AbsTol', 2.24*1e-16);
% 
%     %Flow 1
%     [~, xxM1] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T1],...
%             [X1; reshape(eye(6), [], 1)],...
%             options);
%     Xf1 = (xxM1(end, 1:6))';
%     rhs_f1 = CRTBP(T1, Xf1, mu); 
%     STMf1 = (reshape(xxM1(end, 7:42), 6, 6))';  
% 
%     %Flow 2
%     [~, xxM2] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T2],...
%             [X2; reshape(eye(6), [], 1)],...
%             options);
%     Xf2 = (xxM2(end, 1:6))';
%     rhs_f2 = CRTBP(T2, Xf2, mu); 
%     STMf2 = (reshape(xxM2(end, 7:42), 6, 6))';   %From equations to STM
% 
%     %Flow 3
%     [~, xxM3] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T3],...
%             [X3; reshape(eye(6), [], 1)],...
%             options);
%     Xf3 = (xxM3(end, 1:6))';
%     rhs_f3 = CRTBP(T3, Xf3, mu); 
%     STMf3 = (reshape(xxM3(end, 7:42), 6, 6))';   %From equations to STM
% 
%     %Flow 4
%     [~, xxM4] = ode113(@(t, xM) STMexact_CRTBP(t, xM, mu), [0 T4],...
%             [X4; reshape(eye(6), [], 1)],...
%             options);
%     Xf4 = (xxM4(end, 1:6))';
%     rhs_f4 = CRTBP(T4, Xf4, mu); 
%     STMf4 = (reshape(xxM4(end, 7:42), 6, 6))';   %From equations to STM
% 
%     % Constraint Vector (i.e., positions on nodes)
%     F = [X1(1:3) - X0_des(1:3); ...
%         Xf1(1:3) - X2(1:3); ... 
%         Xf2(1:3) - X3(1:3); ... 
%         Xf3(1:3) - X4(1:3); ... 
%         Xf4(1:3) - Xf_des(1:3)];
% 
%     if nargout > 1
%         %Composition of DF
%         DFrow1 = [eye(3), zeros(3,3), zeros(3,1), zeros(3,6), zeros(3,1),...
%             zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1)];
%         DFrow2 = [STMf1(1:3, :), rhs_f1(1:3), -eye(3), zeros(3,3), zeros(3,1),...
%             zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1)];
%         DFrow3 = [zeros(3,6), zeros(3,1), STMf2(1:3, :), rhs_f2(1:3),...
%             -eye(3), zeros(3,3), zeros(3,1), zeros(3,6), zeros(3,1)];
%         DFrow4 = [zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1),...
%             STMf3(1:3, :), rhs_f3(1:3), -eye(3), zeros(3,3), zeros(3,1)];
%         DFrow5 = [zeros(3,6), zeros(3,1), zeros(3,6), zeros(3,1),...
%             zeros(3,6), zeros(3,1), STMf4(1:3, :), rhs_f4(1:3)];
% 
%         DF = [DFrow1; DFrow2; DFrow3; DFrow4; DFrow5];
%     end
% 
% end



