%===============================================================================%
%   Thesis Title:    Error estimates of a spectral Petrov–Galerkin method for   %
%	                    two-sided fractional reaction–diffusion equations       %
%===============================================================================%
%       PETROV-GALERKIN SPECTRAL METHOD FOR FRACTIONAL REACTION-DIFFUSION       %
%                    WITH JACOBI POLYNOMIALS AND GAUSS-JACOBI QUADRATURE        %
%                                                                               %
%  Author:  Hamidreza Karimi                                                    %
%  Date:    October 2025                                                        %
%  Purpose: Implementation of the Petrov-Galerkin method for:                   %
%           \mathcal{L}_\theta^\alpha u + \mu u = f(x),  x \in (-1,1)           %
%                              f(x)=sin(x)                                      %
%           using Jacobi bases and Gauss-Jacobi quadrature.                     %
%                                                                               %
%  Requires: jags.m, japoly.m, solveEquations.m                                 %
%===============================================================================%

clc; clear; close all;

%% ====================== USER-DEFINED PARAMETERS ======================
mu = 1;                            % Reaction coefficient
N_ref = 512;                       % Reference solution (high-resolution)
N_values = [16, 32, 64, 128];      % Polynomial degrees for convergence study

% CHANGE THIS VALUE TO SWITCH BETWEEN θ = 0.7, 1.0, OR ANY OTHER
theta = 0.7;                       % <--- MODIFY HERE
alpha_list = [1.2, 1.4, 1.6, 1.8]; % α values

%% ====================== AUTO-COMPUTE σ AND σ* ======================
sigma_list = zeros(size(alpha_list));
sigma_star_list = zeros(size(alpha_list));

fprintf('Computing σ and σ* for θ = %.2f using solveEquations.m:\n', theta);
for i = 1:length(alpha_list)
    [sigma_list(i), sigma_star_list(i)] = solveEquations(theta, alpha_list(i));
end

% Plot colors (consistent with thesis figures)
colors = [
    0.00  0.4470  0.7410;   % Blue
    0.8500 0.3250  0.0980;   % Orange
    0.9290  0.6940  0.1250;  % Yellow
    0.4940  0.1840  0.5560   % Purple
];

%% ====================== STORAGE FOR ERRORS & RATES ======================
E1_matrix = zeros(length(alpha_list), length(N_values));
rates     = zeros(length(alpha_list), length(N_values)-1);

%% ====================== MAIN COMPUTATIONAL LOOP ======================
for config = 1:length(alpha_list)
    sigma      = sigma_list(config);
    sigma_star = sigma_star_list(config);
    alpha      = alpha_list(config);
    
    fprintf('\n=== Computing for α = %.1f, θ = %.2f, σ = %.6f, σ* = %.6f ===\n', ...
            alpha, theta, sigma, sigma_star);
    
    % Step 1: Reference solution
    [u_hat_ref, ~, ~] = solve_petrov_galerkin(alpha, sigma, sigma_star, N_ref, @sin, mu);
    
    % Step 2: Compute E1(N)
    for idx = 1:length(N_values)
        N = N_values(idx);
        [u_hat_N, ~, ~] = solve_petrov_galerkin(alpha, sigma, sigma_star, N, @sin, mu);
        
        d_n = u_hat_ref(1:N+1) - u_hat_N;
        d_n = [d_n; u_hat_ref(N+2:end)];
        
        [x_gauss, w_gauss] = jags(2*N_ref + 1, 0, 0);
        u_diff = zeros(size(x_gauss));
        for n = 0:N_ref
            P_n = japoly(n, sigma, sigma_star, x_gauss);
            u_diff = u_diff + d_n(n+1) * P_n;
        end
        
        E1_matrix(config, idx) = sqrt(sum(w_gauss .* (u_diff .^ 2)));
        fprintf('  N = %3d → E1(N) = %.3e\n', N, E1_matrix(config, idx));
    end
    
    % Step 3: Convergence rates
    for i = 1:length(N_values)-1
        rates(config, i) = log2(E1_matrix(config, i) / E1_matrix(config, i+1));
        fprintf('  Rate(N=%d → %d) = %.2f\n', N_values(i), N_values(i+1), rates(config, i));
    end
end

%% ====================== DISPLAY RESULTS TABLE ======================
fprintf('\n\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('                CONVERGENCE RESULTS (θ = %.2f)\n', theta);
fprintf('%s\n', repmat('=', 1, 80));
fprintf('%8s %12s %12s %12s %12s %12s\n', 'α', 'N=16', 'N=32', 'N=64', 'N=128', 'Rate');
fprintf('%s\n', repmat('-', 1, 80));

for config = 1:length(alpha_list)
    fprintf('α=%.1f', alpha_list(config));
    for i = 1:length(N_values)
        fprintf('  %10.2e', E1_matrix(config, i));
    end
    fprintf('  ');
    for i = 1:length(N_values)-1
        fprintf('%.2f ', rates(config, i));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('=', 1, 80));

%% ====================== PLOT ERROR DECAY (LOG-LOG) ======================
figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 12 12]);
hold on; grid on;

for config = 1:length(alpha_list)
    loglog(N_values, E1_matrix(config, :), '-o', ...
        'LineWidth', 2.0, ...
        'Color', colors(config, :), ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', colors(config, :), ...
        'DisplayName', sprintf('\\alpha = %.1f', alpha_list(config)));
end

xlabel('Polynomial Degree N', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('E_1(N)', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Error Decay of Petrov-Galerkin Method (f(x) = sin(x), \\theta = %.2f)', theta), ...
      'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');

set(gca, 'FontSize', 12, 'XScale', 'log', 'YScale', 'log', ...
         'GridAlpha', 0.3, 'Box', 'on');
xticks(N_values);
xticklabels(arrayfun(@num2str, N_values, 'UniformOutput', false));
legend('Location', 'southwest', 'FontSize', 11, 'Interpreter', 'tex');

axis square;
axis tight;


%% ====================== SAVE FIGURE (DYNAMIC FILENAME) ======================
filename = sprintf('Convergence_Theta%.2f.pdf', theta);
set(gcf, 'Toolbar', 'none');
print(gcf, filename, '-dpdf', '-bestfit');
fprintf('\nFigure saved as: %s\n', filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                            AUXILIARY FUNCTIONS                               %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u_hat, x_quad, w_quad] = solve_petrov_galerkin(alpha, sigma, sigma_star, N, f_func, mu)
    k = 0:N;
    const = -sin(pi*alpha) / (sin(pi*sigma) + sin(pi*sigma_star));
    lambda_k = const * exp(gammaln(alpha + k + 1) - gammaln(k + 1));
    
    h_k = exp(log(2)*(sigma + sigma_star + 1) ...
            - log(2*k + sigma + sigma_star + 1) ...
            + gammaln(k + sigma_star + 1) + gammaln(k + sigma + 1) ...
            - gammaln(k + sigma + sigma_star + 1) - gammaln(k + 1));
    
    Q = diag(lambda_k .* h_k);
    
    [x, w] = jags(N+1, alpha, alpha);
    M = zeros(N+1);
    for i = 0:N
        Pk = japoly(i, sigma_star, sigma, x);
        for j = 0:N
            Pn = japoly(j, sigma, sigma_star, x);
            M(i+1,j+1) = sum(w .* Pn .* Pk);
        end
    end
    
    [xx, ww] = jags(N+1, sigma_star, sigma);
    f_vec = zeros(N+1,1);
    for k = 0:N
        Pk = japoly(k, sigma_star, sigma, xx);
        f_vec(k+1) = sum(f_func(xx) .* Pk .* ww);
    end
    
    A = Q + mu * M;
    u_hat = A \ f_vec;
    
    x_quad = xx; w_quad = ww;
end