% Plot all noise covariance matrices (e.g., original, MVGC LWR, MVGC OLS, HMM full, HMM diag) in a single plot. 
% Also compute R-Squared values computed between the original noise covariance matrix and the estimated matrices from the models.
% 
% Inputs:
%   all_covs   - (1, n_models) cell array containing (n_chan, n_chan) noise covariance matrices estimated from the 5 different models
%   titles     - array of strings

function plt_all_noisecovs(all_covs, titles)
    n_states = length(all_covs{1});
    n_models = length(all_covs) - 1; % first in the list is the original 
    r_squareds = cell(n_states, n_models); % size: n_states, n_models
    for istate = 1:n_states
        for imods = 1:n_models
            r_squareds{istate, imods} = r_squared(all_covs{1}{istate}(:), all_covs{imods+1}{istate}(:));
        end
    end

    figure;
    t = tiledlayout(n_states, length(titles));
    title(t, 'Noise covariance matrices')
    for istate = 1:n_states
        for iplot = 1:length(titles)
            nexttile
            imagesc(all_covs{iplot}{istate}); colorbar;
            if iplot == 1 % if it's 'original'
                title(titles{iplot})
            else
                title(strcat(titles{iplot}, sprintf(", r-squared = %.4g", r_squareds{istate, iplot-1})))
            end
        end
    end
end