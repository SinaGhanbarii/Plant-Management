clear, clc
% Define the experimental data
Dp = [1.5, 4, 7.5, 12.5, 20, 27.5, 32.5, 37.5];
freq = [2.07, 2.17, 12.28, 17.23, 37.00, 11.84, 5.44, 11.97];

% Normalize the frequency data
freq_norm = freq / sum(freq);

% Set initial guess for parameters
initial_guess = [15.5, 9];  % Initial guess for [mu, beta]

% Perform optimization
options = optimset('Display', 'iter');
[optimal_params, fval] = fminsearch(@error_func, initial_guess, options);

% Extract optimal parameters
mu_optimal = optimal_params(1);
beta_optimal = optimal_params(2);

% Print results
fprintf('Optimal mu: %.4f\n', mu_optimal);
fprintf('Optimal beta: %.4f\n', beta_optimal);

% Plot the results
x_fine = linspace(0, 40, 1000);
y_fitted = model_func(optimal_params, x_fine);

figure;
bar(Dp, freq_norm);
hold on;
plot(x_fine, y_fitted, 'r-', 'LineWidth', 2);
xlabel('Particle Size (μm)');
ylabel('Normalized Frequency');
legend('Experimental Data', 'Fitted Model');
title('Particle-Size Distribution: Data vs Fitted Model');

% Define the model function
function y = model_func(params, x)
    mu = params(1);
    beta = params(2);
    z = (x - mu) / beta;
    y = (1/beta) * exp(-(z + exp(-z)));
end

% Define the error function to minimize
function err = error_func(params)
    global Dp freq_norm
    predicted = model_func(params, Dp);
    err = sum((predicted - freq_norm).^2);
end