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
%                              f(x) = |sin(x)|                                  %
%           using Jacobi bases and Gauss-Jacobi quadrature.                     %
%                                                                               %
%  Requires: jags.m, japoly.m, solveEquations.m                                 %
%===============================================================================%

clc; clear; close all;

%% ====================== USER-DEFINED PARAMETERS ======================
N_ref = 512;
N_values = [16; 32; 64; 128];

% CHANGE HERE
theta = 0.7;                       % <--- MODIFY θ HERE
alpha_list = [1.2, 1.4, 1.6, 1.8];  % α values

%% ====================== AUTO-COMPUTE σ AND σ* ======================
fprintf('Computing σ and σ* for θ = %.2f...\n', theta);
sigma_list = zeros(size(alpha_list));
sigma_star_list = zeros(size(alpha_list));

for i = 1:length(alpha_list)
    [sigma_list(i), sigma_star_list(i)] = solveEquations(theta, alpha_list(i));
end

fprintf('Done. σ and σ* computed:\n');
for i = 1:length(alpha_list)
    fprintf('  α = %.1f → σ = %.15f, σ* = %.15f\n', alpha_list(i), sigma_list(i), sigma_star_list(i));
end

%% ====================== STORAGE & COLORS ======================
num_cases = length(sigma_list);
E1_values_all = zeros(length(N_values), num_cases);
Rate_all = zeros(length(N_values)-1, num_cases);
colors = [0 0 1; 1 0 0; 0 1 0; 1 0 1];  % Blue, Red, Green, Purple

%% ====================== MAIN LOOP ======================
for case_idx = 1:num_cases
    sigma = sigma_list(case_idx);
    sigma_star = sigma_star_list(case_idx);
    alpha = sigma + sigma_star;
    u_hat_values = cell(length(N_values)+1, 1);

    fprintf('\n--- Case %d: α = %.1f, σ = %.6f, σ* = %.6f ---\n', case_idx, alpha, sigma, sigma_star);

    % Solve for N_ref and all N_values
    for idx = 0:length(N_values)
        if idx == 0
            N = N_ref;
        else
            N = N_values(idx);
        end

        k = 0:N;

        % λ_θ,k^α
        const = -sin(pi*alpha) / (sin(pi*sigma) + sin(pi*sigma_star));
        lambda_k = const * exp(gammaln(alpha + k + 1) - gammaln(k + 1));

        % h_k
        h_k = exp(log(2)*(sigma + sigma_star + 1) ...
                - log(2*k + sigma + sigma_star + 1) ...
                + gammaln(k + sigma + 1) + gammaln(k + sigma_star + 1) ...
                - gammaln(k + sigma + sigma_star + 1) - gammaln(k + 1));

        Q = diag(lambda_k .* h_k);

        % M matrix
        [x, w] = jags(N+1, alpha, alpha);
        M = zeros(N+1);
        for i = 0:N
            Pk = japoly(i, sigma_star, sigma, x);
            for j = 0:N
                Pn = japoly(j, sigma, sigma_star, x);
                M(i+1,j+1) = sum(w .* Pn .* Pk);
            end
        end

        % RHS
        [xx, ww] = jags(N+1, sigma_star, sigma);
        R = zeros(N+1,1);
        for k = 0:N
            Pk = japoly(k, sigma_star, sigma, xx);
            R(k+1) = sum(abs(sin(xx)) .* Pk .* ww);
        end

        % Solve
        u_hat_values{idx+1} = (M + Q) \ R;
    end

    % --- E1(N): FAST & 100% ACCURATE ---
    u_hat_ref = u_hat_values{1};

    for idx = 1:length(N_values)
        N = N_values(idx);
        [x_gauss, w_gauss] = jags(N_ref + 1, 0, 0);
        w_gauss = w_gauss(:).';  % Row vector
        x_gauss = x_gauss(:).';

        % Precompute P_n
        P_all = zeros(N_ref + 1, length(x_gauss));
        for n = 0:N_ref
            P_all(n+1, :) = japoly(n, sigma, sigma_star, x_gauss);
        end

        k = N_ref - N;
        I_sub = zeros(k, k);

        for ii = 1:k
            n_idx = N + ii;
            Pn = P_all(n_idx + 1, :);
            for jj = ii:k
                m_idx = N + jj;
                Pm = P_all(m_idx + 1, :);
                val = w_gauss * (Pn .* Pm).';  % Scalar
                I_sub(ii, jj) = val;
                I_sub(jj, ii) = val;
            end
        end

        tail = u_hat_ref(N + 2:end);
        E1 = sqrt(tail.' * I_sub * tail);
        E1_values_all(idx, case_idx) = E1;

        fprintf('  N = %3d → E1(N) = %.3e\n', N, E1);
    end

    % Convergence rate
    Rate = log(E1_values_all(1:end-1, case_idx) ./ E1_values_all(2:end, case_idx)) / log(2);
    Rate_all(:, case_idx) = Rate;
end

%% ====================== PLOT 1: LOG-LOG ERROR ======================
figure('Color','w','Units','centimeters','Position',[5 5 12 12]);
hold on;
for c = 1:num_cases
    loglog(N_values, E1_values_all(:,c), '-o', ...
        'LineWidth',2, ...
        'Color', colors(c,:), ...
        'MarkerSize',8, ...
        'MarkerFaceColor', colors(c,:));
end
grid on; box on;
set(gca, 'XScale','log','YScale','log','FontSize',12);
xlabel('N','FontSize',14,'FontWeight','bold');
ylabel('E_1(N)','FontSize',14,'FontWeight','bold');
title('Error Decay for f(x) = |sin(x)|','FontSize',15,'FontWeight','bold');
xticks(N_values); xticklabels({'16','32','64','128'});
legend(arrayfun(@(s,ss) sprintf('\\sigma=%.2f, \\sigma^*=%.2f', s, ss), ...
    sigma_list, sigma_star_list, 'UniformOutput', false), ...
    'Location','southwest', 'Interpreter','tex', 'FontSize',11);
axis square;

%% ====================== PLOT 2: CONVERGENCE RATE ======================
figure('Color','w','Units','centimeters','Position',[5 5 12 12]);
hold on;
markers = {'o','s','^','d'};
for c = 1:num_cases
    plot(N_values(2:end), Rate_all(:,c), ['-', markers{c}], ...
         'Color', colors(c,:), ...
         'LineWidth',1.5, ...
         'MarkerSize',8, ...
         'MarkerFaceColor', colors(c,:), ...
         'DisplayName', sprintf('$(\\sigma,\\sigma^*)=(%.3f,%.3f)$', sigma_list(c), sigma_star_list(c)));
end
grid on; box on;
xlabel('$N$','Interpreter','latex','FontSize',14,'FontWeight','bold');
ylabel('$p_{E_1}$','Interpreter','latex','FontSize',14,'FontWeight','bold');
title('Convergence Rate $p_{E_1}$','Interpreter','latex','FontSize',15,'FontWeight','bold');
legend('Location','northwest','Interpreter','latex','FontSize',11);
set(gca,'XScale','log','FontSize',12,'TickLabelInterpreter','latex');
xticks(N_values(2:end)); xticklabels({'32','64','128'});
axis square;
hold off;

%% ====================== FINAL RESULTS TABLE ======================
fprintf('\n\n');
fprintf('%s\n', repmat('=', 1, 90));
fprintf('                     FINAL CONVERGENCE RESULTS — f(x) = |sin(x)| (θ = %.2f)\n', theta);
fprintf('%s\n', repmat('=', 1, 90));
fprintf('%8s %12s %12s %12s %12s %12s %12s\n', ...
        'α', 'N=16', 'N=32', 'N=64', 'N=128', 'Rate16→32', 'Rate32→64', 'Rate64→128');
fprintf('%s\n', repmat('-', 1, 90));

for config = 1:length(alpha_list)
    fprintf('α=%.1f', alpha_list(config));
    for i = 1:length(N_values)
        fprintf('  %10.2e', E1_values_all(i, config));
    end
    fprintf('  ');
    for i = 1:length(N_values)-1
        fprintf('%.2f ', Rate_all(i, config));
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('=', 1, 90));

fprintf('\nAll done! Plots and table generated.\n');