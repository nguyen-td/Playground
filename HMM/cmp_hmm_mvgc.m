% Compare MVAR estimation of 1-state HMM with tsdata_to_var, parameter values are the same as in train_hmmmar.m
%% Add paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/plugins/roiconnect/libs/mvgc_v1.0/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 3/HMM-MAR/'))

%% Load ground truth AR data
load(strcat('outputs/gen_110_eeg.mat'))
X = out_genar{1}';
Arsig = out_genar{2};
x = out_genar{3};
covx = cov(x'); % noise covariance matrix of ground truth data
lambdamax = out_genar{4};

%% Get HMM noise covariance matrices and regression coefficients
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM'

covtypes = {'full', 'diag'};
Ws = cell(1, 2);
covmats = cell(1, 2);

for icov = 1:length(covtypes)
    load(strcat('outputs/hmm_110_', covtypes{icov}, '_eeg.mat'))
    load(strcat('outputs/Gamma_110_', covtypes{icov}, '_eeg.mat'))
    k = size(Gamma, 2); % number of states
    W = getMARmodel(hmm, 1); % get regression coefficient
    Ws{icov} = W;
    
    if strcmpi(covtypes{icov}, 'full')
        order = size(hmm.state(1).Omega.Gam_irate, 1);
        [covmat, ~, ~, ~] = getFuncConn(hmm, k); % noise covariance matrix
        covmats{icov} = covmat;
    else % strcmpi(icov, 'diag')
        covmat = diag(hmm.state(k).Omega.Gam_rate) / (hmm.state(k).Omega.Gam_shape-1); % noise covariance matrix
        covmats{icov} = covmat;
    end
end

%% Fit MVGC MVAR models
% reshape time series back to being multi-trial
n_chans = size(X, 1);
n_trials = 30;
n_samples = size(X, 2) / n_trials;
X = reshape(X', n_chans, n_samples, n_trials);

p = order;
regmodes = {'OLS', 'LWR'};
[A_ols, SIG_ols, E_ols] = tsdata_to_var(X, p, regmodes{1});
[A_lwr, SIG_lwr, E_lwr] = tsdata_to_var(X, p, regmodes{2});

%% Rename noise covariance matrices for clarity
cov_orig = covx;
cov_mvgc_lwr = SIG_lwr;
cov_mvgc_ols = SIG_ols;
cov_hmmfull = covmats{1};
cov_hmmdiag = covmats{2};

all_covs = {cov_orig, cov_mvgc_lwr, cov_mvgc_ols, cov_hmmfull, cov_hmmdiag};

%% Compute correlation coefficients
r_squareds = cell(1, length(all_covs)-1);
for icorr = 1:length(r_squareds)
    r_squareds{icorr} = r_squared(cov_orig(:), all_covs{icorr+1}(:));
end

%% Plot covariance matrices
titles = ["Original MAR model", "MVGC (LWR)", "MVGC (OLS)", "HMM (full)", "HMM (diag)"];
figure;
t = tiledlayout(1, 5);
title(t, 'Noise covariance matrices')
for iplot = 1:length(titles)
    nexttile
    if iplot > 1
        imagesc(all_covs{iplot}); colorbar;
        title(strcat(titles{iplot}, sprintf(", r-squared = %.4g", corrcoefs{iplot-1})))
    else
        imagesc(all_covs{iplot}); colorbar;
        title(titles{iplot})
    end
end

%% Rename regresssion coefficient matrices for clarity
W_orig = Arsig;
W_mvgc_lwr = A_lwr;
W_mvgc_ols = A_ols;
W_hmmfull = reshape(Ws{1}', [], order, order);
W_hmmdiag = reshape(Ws{2}', [], order, order);

all_Ws = {W_orig, W_mvgc_lwr, W_mvgc_ols, W_hmmfull, W_hmmdiag};
min_val = min(cell2mat(all_Ws), [], 'all');
max_val = max(cell2mat(all_Ws), [], 'all');

%% Plot regression coefficient matrices
for iplot = 1:length(titles)
    figure;
    t = tiledlayout(2, 5);
    title(t, titles{iplot})
    for iorder = 1:order
        nexttile
        imagesc(all_Ws{iplot}(:, :, iorder)); colorbar; clim([min_val max_val]);
    end
end

%% Compute PSD of time series
% Code taken from plt_simdata, expanded to also compare original PSD to MVGC PSD