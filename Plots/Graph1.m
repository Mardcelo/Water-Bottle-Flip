% Define the function for hCM/H
hCM_over_H = @(M, f) (1/2) * ((1 + M * f.^2) ./ (1 + M * f));

% Set the value of M and a fine range of f values for higher accuracy
M = 20;
f_values = linspace(0, 1, 1000);  % Increase number of points to 1000 for better resolution

% Calculate hCM/H for the given M and f values
hCM_over_H_values = hCM_over_H(M, f_values);

% Find the minimum value of hCM/H and its corresponding f value
[min_hCM_over_H, min_f_index] = min(hCM_over_H_values);
min_f = f_values(min_f_index);

% Plot the function hCM/H versus f
figure;
plot(f_values, hCM_over_H_values, 'r-', 'LineWidth', 2);
hold on;

% Mark the minimum point on the plot
plot(min_f, min_hCM_over_H, 'ro', 'MarkerSize', 10);

% Label the axes
xlabel('f (filling fraction)');
ylabel('h_{CM}/H', 'Color', 'r');

% Set y-axis limits based on the image provided
ylim([0.1, 0.5]);  % Adjust if needed based on your data



% Add a title
title('h_{CM}/H when M = 20');

% Display the minimum value with high precision
disp(['Minimum h_{CM}/H: ', num2str(min_hCM_over_H, 10)]);
disp(['Corresponding f value: ', num2str(min_f, 10)]);
