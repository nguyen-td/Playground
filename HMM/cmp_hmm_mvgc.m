% Compare MVAR estimation of 1-state HMM with tsdata_to_var, parameter values are the same as in train_hmmmar.m
%% Add paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/plugins/roiconnect/libs/mvgc_v1.0/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 3/HMM-MAR/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'))
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM'


%% Load ground truth AR data
f_path_gtar = 'outputs/gen_110_eeg.mat';
[X_orig, W_orig, noise_orig, ~] = load_groundtruth_ar(f_path_gtar);
cov_orig = cov(noise_orig');

%% Get HMM noise covariance matrices and regression coefficients
n_states = 1; 
n_lags = 10;
[Ws, covmats] = get_W_Sigma_hmm(n_states, n_lags);

%% Fit MVGC MVAR models
n_chans = size(X_orig, 1);
n_trials = 30;
n_samples = size(X_orig, 2) / n_trials;
X_orig = reshape(X_orig, n_chans, n_samples, n_trials); % reshape time series back to being multi-trial

regmodes = {'OLS', 'LWR'};
[W_mvgc_ols, cov_mvgc_ols, ~] = tsdata_to_var(X_orig, n_lags, regmodes{1});
[W_mvgc_lwr, cov_mvgc_lwr, ~] = tsdata_to_var(X_orig, n_lags, regmodes{2});

%% Compute correlation coefficients
cov_hmmfull = covmats{1};
cov_hmmdiag = covmats{2};
all_covs = {cov_orig, cov_mvgc_lwr, cov_mvgc_ols, cov_hmmfull, cov_hmmdiag};

r_squareds = cell(1, length(all_covs)-1);
for icorr = 1:length(r_squareds)
    r_squareds{icorr} = r_squared(cov_orig(:), all_covs{icorr+1}(:));
end

%% Plot covariance matrices
titles = ["Original MAR model", "MVGC (LWR)", "MVGC (OLS)", "HMM (full)", "HMM (diag)"];
plt_all_noisecovs(all_covs, r_squareds, titles)

%% Plot regression coefficient matrices
W_hmmfull = reshape(Ws{1}', n_chans, n_chans, n_lags);
W_hmmdiag = reshape(Ws{2}', n_chans, n_chans, n_lags);
all_Ws = {W_orig, W_mvgc_lwr, W_mvgc_ols, W_hmmfull, W_hmmdiag};

plt_all_Ws(all_Ws, n_lags, titles)
%% Compute PSD of time series, based on plt_simdata.m
Xs_sim = cell(1, 5); % store all AR-generated (simulated) time series

% generate simulated data using MVGC toolbox
X_lwr = reshape(var_to_tsdata(W_mvgc_lwr, cov_mvgc_lwr, n_samples, n_trials), n_chans, []); % remark: generates multi-trial data, gen_2ar.m does not
X_ols = reshape(var_to_tsdata(W_mvgc_ols, cov_mvgc_ols, n_samples, n_trials), n_chans, []);

% generate simulated data using estimated HMM-MAR model
covtypes = {'full', 'diag'};
load outputs/T
for icov = 1:length(covtypes)
    load(strcat('outputs/hmm_110_', covtypes{icov}, '_eeg.mat'))
    X_sim = simhmmmar(T, hmm);
    Xs_sim{icov+3} = X_sim';
end
% Xs_sim order: original, MVGC LWR, MVGC OLS, HMM full, HMM diag
Xs_sim{1} = reshape(X_orig, n_chans, []);
Xs_sim{2} = X_lwr;
Xs_sim{3} = X_ols;

X_pxxs = cell(1, 5);
for ipxx = 1:length(X_pxxs)
    X_pxxs{ipxx} = pwelch(Xs_sim{ipxx}', 100, 50);
    X_pxxs{ipxx} = X_pxxs{ipxx}(2:end-1, :);
end

%% Plot PSDs 
plt_all_psds(X_pxxs, n_chans, titles)