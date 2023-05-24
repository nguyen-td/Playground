%% Load model
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM'
covtype = 'full';
load(strcat('outputs/hmm_110_', covtype, '_eeg.mat'))
load(strcat('outputs/gamma_110_', covtype, '_eeg.mat'))
% load outputs/hmm_33_full_fmri.mat
% load outputs/gamma_33_full_fmri.mat
data_pnts = 500; % plot only the first 500 time points

%% Generate new data
load outputs/T
[X_sim, ~, ~] = simhmmmar(T, hmm);
load(strcat('outputs/gen_110_', covtype, '_eeg.mat'))
X = out_genar{1};
covtype = hmm.train.covtype; % adjust

figure;
tiledlayout(3, 1)

nexttile
plot(X(1:data_pnts, 1))
hold on;
plot(X_sim(1:data_pnts, 1))
legend('original','simulated')
title(strcat("Original and HMM-simulated data, covtype='", covtype,"'"))

nexttile
plot(X(1:data_pnts, 1))
title(strcat("Original data, covtype='", covtype,"'"))

nexttile
plot(X_sim(1:data_pnts, 1), 'color', '#D95319')
title(strcat("HMM-simulated data, covtype='", covtype,"'"))

%% Compute autocorrelation of time series
[X_autocorrx, X_lags] = xcorr(X(1:data_pnts, 1));
[Xsim_autocorr, Xsim_lags] = xcorr(X_sim(1:data_pnts, 1));

figure;
t = tiledlayout(3, 1);
xlabel(t, 'Time lag')
title(t, 'Autocorrelation')

nexttile
plot(X_lags, X_autocorrx)
hold on;
plot(Xsim_lags, Xsim_autocorr, 'color', '#D95319')
title(strcat("Original and HMM-simulated data, covtype='", covtype,"'"))
legend('original', 'simulated')

nexttile
plot(X_lags, X_autocorrx)
title(strcat("Original data, covtype='", covtype,"'"))

nexttile
plot(Xsim_lags, Xsim_autocorr, 'color', '#D95319')
title(strcat("HMM-simulated data, covtype='", covtype, "'"))

%% Compute PSD of time series
X_pxx = pwelch(X, 100, 50);
Xsim_pxx = pwelch(X_sim, 100, 50);

% remove first and last frequencies due to edge artifacts
X_pxx = X_pxx(2:end-1, :);
Xsim_pxx = Xsim_pxx(2:end-1, :);

figure;
t = tiledlayout(2, 5);
title(t, strcat("PSD of original and HMM-simulated data, covtype='", covtype, "'"))

for i = 1:size(X_pxx, 2)
    x = 1:size(X_pxx, 1);
    R = corrcoef(log(abs(X_pxx(:, i))), log(abs(Xsim_pxx(:, i))));

    nexttile
    plot(x+1, log(abs(X_pxx(:, i))))
    hold on;
    plot(x+1, log(abs(Xsim_pxx(:, i))), 'color', '#D95319')
    ylabel('log(abs(Power))')
    legend('original', 'simulated')
    title(['r = ' num2str(R(1, 2))]) % pick any off-diagonal of the matrix, cf. https://de.mathworks.com/help/matlab/ref/corrcoef.html#bunkaln
end