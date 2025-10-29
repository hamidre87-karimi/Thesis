%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                  demo_sigma_computation.m                                    %
%         Generate Table of σ and σ* for Various θ and α                       %
%                                                                              %
%  Author:  Hamidreza Karimi                                                   %
%  Date:    October 2025                                                       %
%                                                                              %
%  Purpose: Demonstrate the computation of Petrov–Galerkin parameters        %
%            σ and σ* for different values of θ and α.                        %
%            Results are displayed and saved as 'sigma_table.csv'.           %
%                                                                              %
%  Requires: solveEquations.m                                                  %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% ====================== INPUT PARAMETERS ======================
theta_values = [0.7, 1.0];           % Test values for θ
alpha_values = [1.2, 1.4, 1.6, 1.8];  % Test values for α

%% ====================== COMPUTE σ AND σ* ======================
results = table();  % Initialize empty table

fprintf('Computing σ and σ* for θ ∈ [%.1f, %.1f] and α ∈ [%.1f, %.1f]...\n', ...
        min(theta_values), max(theta_values), min(alpha_values), max(alpha_values));

for i = 1:length(theta_values)
    for j = 1:length(alpha_values)
        theta = theta_values(i);
        alpha = alpha_values(j);
        
        % Solve nonlinear system
        [sigma, sigma_star] = solveEquations(theta, alpha);
        
        % Append to table
        results = [results; table(theta, alpha, sigma, sigma_star)];
    end
end

%% ====================== DISPLAY RESULTS ======================
disp(' ');
disp('=== σ AND σ* COMPUTATION RESULTS ===');
disp(results);


