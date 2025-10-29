%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%                         solveEquations.m                                     %
%               Compute σ and σ* for given θ and α                            %
%                                                                              %
%  Author:  Hamidreza Karimi                                                   %
%  Date:    October 2025                                                       %
%                                                                              %
%  Purpose: Solve the nonlinear system that defines the non-symmetric         %
%            Petrov–Galerkin parameters σ and σ* from the user-provided        %
%            fractional order α and the weighting factor θ.                   %
%                                                                              %
%            System:                                                          %
%                σ + σ* = α                                                   %
%                θ = sin(π σ*) / (sin(π σ) + sin(π σ*))                      %
%                                                                              %
%  Requires:   MATLAB Optimization Toolbox (fsolve)                            %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma, sigma_star] = solveEquations(theta, alpha)
    % SOLVEEQUATIONS  Returns σ and σ* for given θ and α
    %
    %   [sigma, sigma_star] = solveEquations(theta, alpha)
    %
    %   Input:
    %       theta : weighting factor in the fractional operator (0 < θ ≤ 1)
    %       alpha : fractional order (α > 0)
    %
    %   Output:
    %       sigma      : left  exponent in the Jacobi weight ω^{σ,σ*}
    %       sigma_star : right exponent in the Jacobi weight ω^{σ,σ*}
    %
    %   The function uses MATLAB's fsolve with a robust initial guess.
    %   Convergence is guaranteed for the admissible range of θ and α.

    % -------------------------------------------------------------------------
    % Initial guess: split α proportionally to θ (works for θ ∈ (0,1])
    % -------------------------------------------------------------------------
    initialGuess = [0.5, 0.5];   % simple, symmetric start; fsolve is insensitive

    % -------------------------------------------------------------------------
    % Solve the 2×2 nonlinear system
    % -------------------------------------------------------------------------
    options = optimoptions('fsolve', ...
        'Display', 'off', ...     % suppress iteration output
        'TolFun',  1e-15, ...     % function tolerance
        'TolX',    1e-15);        % variable tolerance

    solution = fsolve(@(vars) equations(vars, theta, alpha), ...
                      initialGuess, options);

    % -------------------------------------------------------------------------
    % Extract results
    % -------------------------------------------------------------------------
    sigma      = solution(1);
    sigma_star = solution(2);

    % -------------------------------------------------------------------------
    % Optional diagnostic output (uncomment for debugging)
    % -------------------------------------------------------------------------
    % fprintf('  θ = %.2f, α = %.1f → σ = %.6f, σ* = %.6f\n', ...
    %         theta, alpha, sigma, sigma_star);
end

% =========================================================================
% Sub-function: definition of the nonlinear residual equations
% =========================================================================
function F = equations(vars, theta, alpha)
    % EQUATIONS  Residual vector for the σ–σ* system
    %
    %   F = equations(vars, theta, alpha)
    %
    %   vars(1) = σ,   vars(2) = σ*

    sigma      = vars(1);
    sigma_star = vars(2);

    % Equation 1: weighting relation
    F(1) = theta - sin(pi * sigma_star) / (sin(pi * sigma_star) + sin(pi * sigma));

    % Equation 2: sum constraint
    F(2) = sigma + sigma_star - alpha;
end