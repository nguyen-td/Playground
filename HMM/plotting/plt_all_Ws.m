% Plot all correlation coefficient matrices (original, MVGC LWR, MVGC OLS, HMM full, HMM diag) in individual plots.
%
% Inputs:
%   all_Ws - (1, 5) cell array containing (n_chan, n_chan, n_lags) regression coefficient matrices
%   n_lags - number of lags
%   titles - array of strings

function plt_all_Ws(all_Ws, n_lags, titles)
    min_val = min(cell2mat(all_Ws), [], 'all');
    max_val = max(cell2mat(all_Ws), [], 'all');
    for iplot = 1:length(titles)
        figure;
        t = tiledlayout(1, 10);
        title(t, titles{iplot})
        for ilags = 1:n_lags
            nexttile
            imagesc(all_Ws{iplot}(:, :, ilags)); colorbar; clim([min_val max_val]);
        end
    end
end