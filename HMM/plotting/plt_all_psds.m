% Plot all PSDs of AR-generated time series (original and all estimated ones) in a single plot. 
% Also compute Pearson's correlation coefficient.
%
% Input:
%   X_pxxs  - (1, 5) cell array containing (n, n_chan) power time series in the following order: original, MVGC LWR, MVGC OLS, HMM full, HMM diag
%   n_chans - number of channels
%   titles  - array of strings

function plt_all_psds(X_pxxs, n_chans, titles)
    r_all = zeros(1, length(X_pxxs) - 1); % store correlation coefficients
    figure;
    t = tiledlayout(2, n_chans/2);
    title(t, "PSDs of original and estimated data")

    for ichans = 1:n_chans
        nexttile
        
        % plot original PSD
        x = 1:size(X_pxxs{1}, 1);
        plot(x+1, log(abs(X_pxxs{1}(:, ichans))), 'LineWidth', 2) 
        ylabel('log(abs(Power))')
        hold on;

        for idata = 2:length(X_pxxs)
            R = corrcoef(log(abs(X_pxxs{1}(:, ichans))), log(abs(X_pxxs{idata}(:, ichans)))); % correlation between original PSD and estimated PSD
            r_all(idata-1) = R(1, 2);

            plot(x+1, log(abs(X_pxxs{idata}(:, ichans))))
        end
%         legend('original', ...
%               ['MVGC LWR, r = ' num2str(r_all(1))], ...
%               ['MVGC OLS, r = ' num2str(r_all(2))], ...
%               ['HMM (full), r = ' num2str(r_all(3))], ...
%               ['HMM (diag), r = ' num2str(r_all(4))])     
        leg_titles = cell(length(titles), 1);
        for ileg = 1:length(leg_titles)
            if ileg == 1
                leg_titles{ileg} = titles(ileg);
            else
                leg_titles{ileg} = strcat(titles(ileg), ", r = ", num2str(r_all(ileg - 1)));
            end
        end
        legend(leg_titles)
    end
end