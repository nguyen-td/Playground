%% Load model
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM'
load outputs/hmm_33_full_fmri.mat
load outputs/gamma_33_full_fmri.mat
data_pnts = 500; % plot only the first 500 time points

%% Generate new data
load outputs/T_fmri
[X_sim, ~, ~] = simhmmmar(T, hmm);
load outputs/X_fmri
covtype = 'full'; % adjust

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
[Xsim_autocorr, Xsim_lags] = xcorr(X(1:data_pnts, 1));

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
title(strcat("HMM-simulated data, covtype='", covtype,"'"))

%% Compute PSD of time series
X_pxx = pwelch(X(1:data_pnts, 1), data_pnts, floor(data_pnts/2));
Xsim_pxx = pwelch(X_sim(1:data_pnts, 1), data_pnts, floor(data_pnts/2));
% X_pxx = fft(X(1:data_pnts, 1));
% Xsim_pxx = fft(X_sim(1:data_pnts, 1));

figure;
t = tiledlayout(2, 1);
title(t, 'PSD')

nexttile
plot(log(abs(X_pxx)))
title(strcat("Original data, covtype='", covtype,"'"))
ylabel('log(abs(Power))')

nexttile
plot(log(abs(Xsim_pxx)), 'color', '#D95319')
title(strcat("HMM-simulated data, covtype='", covtype,"'"))
ylabel('log(abs(Power))')