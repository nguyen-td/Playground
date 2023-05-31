% Plot all noise covariance matrices (original, MVGC LWR, MVGC OLS, HMM full, HMM diag) in a single plot.
% 
% Inputs:
%   all_covs   - (1, 5) cell array containing (n_chan, n_chan) noise covariance matrices estimated from the 5 different models
%   r_squareds - (1, 4) cell array containing R-Squared values computed between the original noise covariance matrix and the 
%                 estimated matrices from the models. Pass 0 if you did not compute R-Squared values.
%   titles     - array of strings

function plt_all_noisecovs(all_covs, r_squareds, titles)
    figure;
    t = tiledlayout(1, length(titles));
    title(t, 'Noise covariance matrices')
    for iplot = 1:length(titles)
        nexttile
        if iplot > 1
            imagesc(all_covs{iplot}); colorbar;
            if r_squareds == 0
                title(titles{iplot})
            else
                title(strcat(titles{iplot}, sprintf(", r-squared = %.4g", r_squareds{iplot-1})))
            end
        else
            imagesc(all_covs{iplot}); colorbar;
            title(titles{iplot})
        end
    end
end