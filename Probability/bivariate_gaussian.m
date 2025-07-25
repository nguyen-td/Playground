% Parameters
mu = [0 0]; % Mean [mu_x, mu_y]
Sigma = [1 -0.9; -0.9 1]; % Covariance matrix

% Create grid
x1 = linspace(-3, 3, 100);
x2 = linspace(-3, 3, 100);
[X1, X2] = meshgrid(x1, x2);
X = [X1(:) X2(:)];

% Evaluate bivariate normal PDF
F = mvnpdf(X, mu, Sigma);
F = reshape(F, length(x2), length(x1));

% Plot contour
figure;
contour(X1, X2, F, 10); % '10' is number of contour levels
xlabel('X_1');
ylabel('X_2');
title('Bivariate Normal Distribution');
axis equal;
colorbar;
