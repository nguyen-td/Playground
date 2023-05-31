% Compute the coefficient of determination (R-Squared).
%
% Inputs:
%   y - (n, 1) ground truth data
%   f - (n, 1) fitted values
%
% Output:
%   r2 - coefficient of determination [0 < r2 < 1]

function r2 = r_squared(y, f)
    y_bar = mean(y);
    ss_res = (y - f)' * (y - f);
    ss_tot = (y - y_bar)' * (y - y_bar);
    r2 = 1 - ss_res / ss_tot;
end