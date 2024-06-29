% Define the experimental data
Dp = [1.5, 4, 7.5, 12.5, 20, 27.5, 32.5, 37.5]; % Midpoints of intervals
percentages = [2.07, 2.17, 12.28, 17.23, 37.00, 11.84, 5.44, 11.97];
weights = percentages / 100;

% Define the probability distribution function
f = @(Dp, mu, beta) (1/beta) * exp(-((Dp - mu)/beta + exp(-(Dp - mu)/beta)));

% Define the objective function to minimize
objective = @(params) sum((weights - f(Dp, params(1), params(2))).^2);

% Use optimization to find best μ and β
initial_guess = [15, 10]; % Starting guess based on the data range
options = optimset('Display', 'iter');
optimal_params = fminsearch(objective, initial_guess, options);

mu_optimal = optimal_params(1);
beta_optimal = optimal_params(2);

% Calculate z for Dp = 15.43
Dp_given = 15.43;
z = (Dp_given - mu_optimal) / beta_optimal;

% Generate points for plotting the fitted function
Dp_plot = linspace(0, 50, 1000);
f_plot = f(Dp_plot, mu_optimal, beta_optimal);

% Plot the results
figure;
bar(Dp, weights);
hold on;
plot(Dp_plot, f_plot, 'r-', 'LineWidth', 2);
xlabel('Particle Size D_p (μm)');
ylabel('Probability Density');
title('Particle-Size Distribution and Fitted Function');
legend('Experimental Data', 'Fitted Function');
grid on;

% Calculate efficiency factor η (assuming it's the CDF up to Dp_given)
eta = integral(@(x) f(x, mu_optimal, beta_optimal), 0, Dp_given);

% Display results
fprintf('Optimized mu: %.4f\n', mu_optimal);
fprintf('Optimized beta: %.4f\n', beta_optimal);
fprintf('z for Dp = 15.43: %.4f\n', z);
fprintf('Efficiency factor η: %.4f\n', eta);

% Plot the CDF
figure;
Dp_cdf = linspace(0, 50, 1000);
cdf_values = zeros(size(Dp_cdf));
for i = 1:length(Dp_cdf)
    cdf_values(i) = integral(@(x) f(x, mu_optimal, beta_optimal), 0, Dp_cdf(i));
end
plot(Dp_cdf, cdf_values, 'b-', 'LineWidth', 2);
hold on;
plot([Dp_given, Dp_given], [0, eta], 'r--', 'LineWidth', 2);
plot([0, Dp_given], [eta, eta], 'r--', 'LineWidth', 2);
xlabel('Particle Size D_p (μm)');
ylabel('Cumulative Probability');
title('Cumulative Distribution Function');
legend('CDF', 'Dp = 15.43 μm');
grid on;