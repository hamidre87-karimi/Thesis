%===============================================================================%
%   Thesis Title:    Error estimates of a spectral Petrov–Galerkin method for   %
%	                    two-sided fractional reaction–diffusion equations       %
%===============================================================================
%       PETROV-GALERKIN SPECTRAL METHOD FOR FRACTIONAL REACTION-DIFFUSION       %
%                    WITH JACOBI POLYNOMIALS AND GAUSS-JACOBI QUADRATURE        %
%                                                                               %
%  Author:  Hamidreza Karimi                                                    %
%  Date:    October 2025                                                        %
%  Purpose: Implementation of the Petrov-Galerkin method for:                   %
%           \mathcal{L}_\theta^\alpha u + \mu u = f(x),  x \in (-1,1)           %
%                       f(x) = (1-x²)^(-0.4) * sin(x)                           %
%           using Jacobi bases and Gauss-Jacobi quadrature.                     %
%                                                                               %
%  Requires: jags.m, japoly.m, solveEquations.m                                 %
%===============================================================================%

clc; clear; close all;

%% ====================== USER-DEFINED PARAMETERS ======================
N_ref    = 512;                     % Reference solution
N_values = [16, 32, 64, 128];       % Polynomial degrees

theta    = 0.7;                     % <--- MODIFY θ HERE
alpha_list = [1.2, 1.4, 1.6, 1.8];   % α values

%% ====================== AUTO-COMPUTE σ AND σ* ======================
fprintf('Computing σ and σ* for θ = %.2f using solveEquations.m:\n', theta);
sigma_list = zeros(size(alpha_list));
sigma_star_list = zeros(size(alpha_list));

for i = 1:length(alpha_list)
    [sigma_list(i), sigma_star_list(i)] = solveEquations(theta, alpha_list(i));
end

%% ====================== STORAGE & COLORS ======================
num_cases = length(alpha_list);
E2_matrix = zeros(num_cases, length(N_values));
rates     = zeros(num_cases, length(N_values)-1);
colors = [0 0 1; 1 0 0; 0 1 0; 1 0 1];  % Blue, Red, Green, Purple

%% ====================== MAIN LOOP  ======================
for config = 1:num_cases
    sigma = sigma_list(config);
    sigma_star = sigma_star_list(config);
    alpha = alpha_list(config);
    
    fprintf('\n=== Computing for α = %.1f, θ = %.2f, σ = %.6f, σ* = %.6f ===\n', ...
            alpha, theta, sigma, sigma_star);
    
    u_hat_values = cell(length(N_values)+1, 1);

    % Solve for N_ref and all N_values
    for idx = 0:length(N_values)
        if idx == 0
            N = N_ref;
        else
            N = N_values(idx);
        end
        k = 0:N;

        % Stiffness: λ_θ,k^α
        const = -sin(pi*alpha)/(sin(pi*sigma) + sin(pi*sigma_star));
        lambda_k = const * exp(gammaln(alpha + k + 1) - gammaln(k + 1));
        h_k = exp(log(2)*(sigma + sigma_star + 1) ...
                - log(2*k + sigma + sigma_star + 1) ...
                + gammaln(k + sigma + 1) + gammaln(k + sigma_star + 1) ...
                - gammaln(k + sigma + sigma_star + 1) - gammaln(k + 1));
        Q = diag(lambda_k .* h_k);

        % Mass matrix with Jacobi weight
        [x, w_base] = jags(N+1, 0, 0);
        w_jacobi = exp(sigma*log(1-x) + sigma_star*log(1+x));
        w = w_base .* w_jacobi;

        M = zeros(N+1);
        for i = 0:N
            Pk = japoly(i, sigma_star, sigma, x);
            for j = 0:N
                Pn = japoly(j, sigma, sigma_star, x);
                M(i+1,j+1) = sum(w .* Pn .* Pk);
            end
        end

        % RHS
        [xx, ww_base] = jags(N+1, 0, 0);
        ww_jacobi = exp(sigma_star*log(1-xx) + sigma*log(1+xx));
        ww = ww_base .* ww_jacobi;
        R = zeros(N+1,1);
        for k = 0:N
            Pk = japoly(k, sigma_star, sigma, xx);
            f = ((1 - xx.^2).^(-0.4)) .* sin(xx);
            R(k+1) = sum(f .* Pk .* ww);
        end

        u_hat_values{idx+1} = (M + Q) \ R;
    end

    % High-resolution reference
    [xxxxx, wwwww_base] = jags(1024, 0, 0);
    wwwww_jacobi = exp(-sigma*log(1-xxxxx) + -sigma_star*log(1+xxxxx));
    wwwww = wwwww_base .* wwwww_jacobi;

    u_ref = zeros(length(xxxxx), 1);
    for n = 0:N_ref
        P_n = japoly(n, sigma, sigma_star, xxxxx);
        weight = exp(sigma*log(1-xxxxx) + sigma_star*log(1+xxxxx));
        u_ref = u_ref + u_hat_values{1}(n+1) * P_n .* weight;
    end

    % Approximate solutions
    u_N = zeros(length(xxxxx), length(N_values));
    for idx = 1:length(N_values)
        N = N_values(idx);
        for n = 0:N
            P_n = japoly(n, sigma, sigma_star, xxxxx);
            weight = exp(sigma*log(1-xxxxx) + sigma_star*log(1+xxxxx));
            u_N(:, idx) = u_N(:, idx) + u_hat_values{idx+1}(n+1) * P_n .* weight;
        end
    end

    % E2 error
    diff_u = u_ref - u_N;
    E2_vals = sqrt(sum((diff_u.^2) .* wwwww, 1))';
    E2_matrix(config, :) = E2_vals;

    % Display E2(N)
    for i = 1:length(N_values)
        fprintf('  N = %3d → E2(N) = %.3e\n', N_values(i), E2_vals(i));
    end

    % Convergence rate
    for i = 1:length(N_values)-1
        rates(config, i) = log2(E2_vals(i) / E2_vals(i+1));
        fprintf('  Rate(N=%d → %d) = %.2f\n', N_values(i), N_values(i+1), rates(config, i));
    end
end

%% ====================== PLOT 1: ERROR DECAY (LOG-LOG) ======================
figure('Color','w','Units','centimeters','Position',[5 5 12 12]);
hold on; grid on;
for config = 1:num_cases
    loglog(N_values, E2_matrix(config,:), '-o', ...
        'LineWidth',2, 'Color',colors(config,:), ...
        'MarkerSize',8, 'MarkerFaceColor',colors(config,:), ...
        'DisplayName', sprintf('\\alpha = %.1f', alpha_list(config)));
end
xlabel('Polynomial Degree N','FontSize',14,'FontWeight','bold');
ylabel('E_2(N)','FontSize',14,'FontWeight','bold');
title(sprintf('Error Decay for f(x) = (1-x^2)^{-0.4} \\sin(x) (\\theta = %.2f)', theta), ...
      'FontSize',13,'FontWeight','bold','Interpreter','tex');
set(gca,'FontSize',12,'XScale','log','YScale','log','GridAlpha',0.3,'Box','on');
xticks(N_values); xticklabels(arrayfun(@num2str, N_values, 'UniformOutput', false));
legend('Location','southwest','FontSize',11,'Interpreter','tex');
axis square; axis tight;

%% ====================== PLOT 2: CONVERGENCE RATE ======================
figure('Color','w','Units','centimeters','Position',[5 5 12 12]);
hold on; grid on;
markers = {'o','s','^','d'};
for config = 1:num_cases
    plot(N_values(2:end), rates(config,:), ['-', markers{config}], ...
        'Color',colors(config,:), 'LineWidth',1.5, 'MarkerSize',8, ...
        'DisplayName', sprintf('$(\\sigma,\\sigma^*)=(%.3f,%.3f)$', sigma_list(config), sigma_star_list(config)));
end
xlabel('$N$','Interpreter','latex','FontSize',14,'FontWeight','bold');
ylabel('$p_{E_2}$','Interpreter','latex','FontSize',14,'FontWeight','bold');
title(sprintf('Convergence Rate $p_{E_2}$ (\\theta = %.2f)', theta), ...
      'Interpreter','latex','FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',11,'Interpreter','latex');
set(gca,'XScale','log','FontSize',12,'TickLabelInterpreter','latex');
xticks(N_values(2:end)); xticklabels({'32','64','128'});
axis square;

%% ====================== SAVE PLOTS ======================
print(1, 'Error_Decay_E2.pdf', '-dpdf', '-bestfit');
print(2, 'Rate_E2.pdf', '-dpdf', '-bestfit');

%% ====================== FINAL RESULTS TABLE ======================
fprintf('\n\n');
fprintf('%s\n', repmat('=', 1, 90));
fprintf('                FINAL CONVERGENCE RESULTS — E_2 ERROR (θ = %.2f)\n', theta);
fprintf('%s\n', repmat('=', 1, 90));
fprintf('%8s %12s %12s %12s %12s %12s %12s\n', ...
        'α', 'N=16', 'N=32', 'N=64', 'N=128', 'Rate16→32', 'Rate32→64', 'Rate64→128');
fprintf('%s\n', repmat('-', 1, 90));

for config = 1:num_cases
    fprintf('α=%.1f', alpha_list(config));
    for i = 1:length(N_values)
        fprintf('  %10.2e', E2_matrix(config,i));
    end
    fprintf('  ');
    for i = 1:length(N_values)-1
        fprintf('%.2f ', rates(config,i));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('=', 1, 90));

fprintf('\nAll done! Plots and table saved.\n');