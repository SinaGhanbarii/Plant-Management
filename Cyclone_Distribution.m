% Exercise 1.1 - Advanced Numerical Methods
% Fitting particle-size distribution data

% Clear workspace and close all figures
clear all;
close all;
clc;

% Define the experimental data
Dp = [1.5, 4, 7.5, 12.5, 20, 27.5, 32.5, 37.5];  % Particle size (μm)
freq = [2.07, 2.17, 12.28, 17.23, 37.00, 11.84, 5.44, 11.97];  % Frequency (%)

% Normalize the frequency data
freq_norm = freq / sum(freq);

% Define the model function
model_func = @(params, x) (1./params(2)) .* exp(-(((x - params(1))./params(2)) + exp(-(x - params(1))./params(2))));

% Define the error function to minimize
error_func = @(params) sum((model_func(params, Dp) - freq_norm).^2);

% Set initial guess and bounds for parameters [mu, beta]
initial_guess = [15, 9];
lb = [0, 0];  % Lower bounds
ub = [40, 40];  % Upper bounds

% Perform optimization using fmincon (constrained optimization)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[optimal_params, fval, exitflag, output] = fmincon(error_func, initial_guess, [], [], [], [], lb, ub, [], options);

% Extract optimal parameters
mu_optimal = optimal_params(1);
beta_optimal = optimal_params(2);

% Print results
fprintf('Optimization Results:\n');
fprintf('Optimal mu: %.4f μm\n', mu_optimal);
fprintf('Optimal beta: %.4f μm\n', beta_optimal);
fprintf('Final error: %.6e\n', fval);
fprintf('Exit flag: %d\n', exitflag);
fprintf('Number of iterations: %d\n', output.iterations);

% Generate fitted curve
x_fine = linspace(0, 40, 1000);
y_fitted = model_func(optimal_params, x_fine);

% Plot the results
figure('Name', 'Particle-Size Distribution: Data vs Fitted Model', 'NumberTitle', 'off');
bar(Dp, freq_norm, 'FaceColor', [0.8 0.8 0.8]);
hold on;
plot(x_fine, y_fitted, 'r-', 'LineWidth', 2);
xlabel('Particle Size D_p (μm)', 'FontSize', 12);
ylabel('Normalized Frequency', 'FontSize', 12);
legend('Experimental Data', 'Fitted Model', 'Location', 'best');
title('Particle-Size Distribution: Data vs Fitted Model', 'FontSize', 14);
grid on;
hold off;

% Calculate R-squared value
SS_tot = sum((freq_norm - mean(freq_norm)).^2);
SS_res = sum((freq_norm - model_func(optimal_params, Dp)).^2);
R_squared = 1 - SS_res / SS_tot;
fprintf('R-squared value: %.4f\n', R_squared);
