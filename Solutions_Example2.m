%===============================================================================%
%   Thesis Title:    Error estimates of a spectral Petrov–Galerkin method for   %
%	                    two-sided fractional reaction–diffusion equations       %
%===============================================================================%
%                                                                               %
%           APPROXIMATE SOLUTION u_N(x) - PETROV-GALERKIN METHOD                %
%                 Jacobi Bases | Dynamic θ | N = 128                            %
%                                                                               %
%  Author:  Hamidreza Karimi                                                    %
%  Date:    October 2025                                                        %
%  Purpose: Plot u_N(x) for different α values using Petrov-Galerkin method     %
%           with σ and σ* computed automatically from θ and α.                  %
%                              f(x)=|sin(x)|                                    %
%                                                                               %
%  Requires: jags.m, japoly.m, solveEquations.m                                 %
%                                                                               %
%===============================================================================%

clc; clear; close all;

%% ====================== USER-DEFINED PARAMETERS ======================
N = 128;                                      % Polynomial degree
x_plot = linspace(-1, 1, 1000)';              % High-resolution plotting grid

theta = 0.7;                                  % MODIFY θ HERE
alpha_list = [1.2, 1.4, 1.6, 1.8];             % α values

%% ====================== AUTO-COMPUTE σ AND σ* ======================
sigma_list = zeros(size(alpha_list));
sigma_star_list = zeros(size(alpha_list));

fprintf('Computing σ and σ* for θ = %.2f using solveEquations.m:\n', theta);
for i = 1:length(alpha_list)
    [sigma_list(i), sigma_star_list(i)] = solveEquations(theta, alpha_list(i));
end

line_styles = {'-', '--', ':', '-.'};
colors = lines(length(alpha_list));

%% ====================== MAIN LOOP: SOLVE & PLOT ======================
figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 16 10]);
hold on; grid on; box on;

for idx = 1:length(alpha_list)
    sigma      = sigma_list(idx);
    sigma_star = sigma_star_list(idx);
    alpha      = alpha_list(idx);
    
    fprintf('Solving for α = %.1f, θ = %.2f, σ = %.6f, %c%c = %.6f\n', ...
            alpha, theta, sigma, char(963), char(42), sigma_star);
    
    const = -sin(pi * alpha) / (sin(pi * sigma) + sin(pi * sigma_star));
    k = 0:N;
    lambda_k = const * exp(gammaln(alpha + k + 1) - gammaln(k + 1));
    
    h_k = exp(log(2) * (sigma + sigma_star + 1) ...
            - log(2 * k + sigma + sigma_star + 1) ...
            + gammaln(k + sigma + 1) + gammaln(k + sigma_star + 1) ...
            - gammaln(k + sigma + sigma_star + 1) - gammaln(k + 1));
    Q = diag(lambda_k .* h_k);
    
    [x_quad, w_quad] = jags(N+1, alpha, alpha);
    M = zeros(N+1);
    for i = 0:N
        Pk = japoly(i, sigma_star, sigma, x_quad);
        for j = 0:N
            Pn = japoly(j, sigma, sigma_star, x_quad);
            M(i+1,j+1) = sum(w_quad .* Pn .* Pk);
        end
    end
    
    [xx, ww] = jags(N+1, sigma_star, sigma);
    R = zeros(N+1,1);
    for kk = 0:N
        Pkk = japoly(kk, sigma_star, sigma, xx);
        R(kk+1) = sum(abs(sin(xx)) .* Pkk .* ww);
    end
    
    u_hat = (Q + M) \ R;
    
    weight = (1 - x_plot).^sigma .* (1 + x_plot).^sigma_star;
    u_N = zeros(size(x_plot));
    for nn = 0:N
        P_n = japoly(nn, sigma, sigma_star, x_plot);
        u_N = u_N + u_hat(nn+1) * weight .* P_n;
    end
    
    plot(x_plot, u_N, ...
        'LineStyle', line_styles{idx}, ...
        'Color', colors(idx,:), ...
        'LineWidth', 2.2, ...
        'DisplayName', sprintf('\\theta=%.2f, \\alpha=%.1f', theta, alpha));
end

%% ====================== FINALIZE FIGURE ======================
xlabel('$x$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
ylabel('$u_N(x)$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');

legend('Location', 'northwest', 'FontSize', 11, 'Interpreter', 'tex');
set(gca, 'FontSize', 12, 'GridAlpha', 0.3, 'LineWidth', 1.1);

xlim([-1.02, 1.02]);
ylim(ylim + [-0.02, 0.02]);
set(gcf, 'Toolbar', 'none');

fprintf('\nPlot displayed. Close the figure to continue.\n');