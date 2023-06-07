% Compare MVAR estimation of 3-state HMM to the simulated state-dependent MVAR model in get_state_ar.m.

%% Add paths and set some parameters
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 3/HMM-MAR/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'))
DIROUT = 'outputs/';
n_states = 3;
P = 10; % number of lags

%% Generate simulated data using estimated HMM-MAR model
load(sprintf(strcat(DIROUT, 'data_sim_gt_%d%d.mat'), n_states, P))

covtypes = {'full', 'diag'};
X = cell(1, length(covtypes) + 1);
X{1} = data_gt; % simulated ground truth data

load outputs/T
for icov = 1:length(covtypes)
    load(strcat('outputs/hmm_310_', covtypes{icov}, '_eeg.mat'))
    X{icov + 1} = simhmmmar(T, hmm);
end
n_chans = size(data_gt, 2);

%% Compute PSD of time series
X_pxxs = cell(1, length(covtypes) + 1);
for ipxx = 1:length(X_pxxs)
    X_pxxs{ipxx} = pwelch(X{ipxx}, 100, 50);
end

%% Plot PSDs 
titles = ["original", "HMM (full)", "HMM (diag)"];
plt_all_psds(X_pxxs, n_chans, titles)

%% Get noise covariance matrices and regression coefficients
load(sprintf(strcat(DIROUT, 'Sigma_sim_gt_%d%d.mat'), n_states, P)); % load ground truth matrices
[Ws, covmats] = get_W_Sigma_hmm(n_states, P); % load estimated matrices

cov_orig = Sigma;
cov_hmmfull = covmats(1, :);
cov_hmmdiag = covmats(2, :);
all_covs = {cov_orig, cov_hmmfull, cov_hmmdiag};

%% Plot noise covariance matrices
plt_all_noisecovs(all_covs, titles)

%% Get Viterbi paths
vpaths = cell(1, length(covtypes));
for icov = 1:length(covtypes)
    load(strcat('outputs/hmm_310_', covtypes{icov}, '_eeg.mat'))
    vpaths{icov} = hmmdecode(X{1}, T, hmm, 1);
end

%% Plot state time series
fs = 200;
for icov  = 1:length(covtypes)
    load(strcat('outputs/Gamma_310_', covtypes{icov}, '_eeg.mat'))
    t = (1:size(Gamma, 1)) / fs;
    figure; plot(vpaths{1})
    title(strcat("Viterbi path, covtype = '", covtypes{icov}, "'")) 

    figure; imagesc(Gamma'); colorbar;
    title(strcat("State time series Gamma, covtype = '", covtypes{icov}, "'")) 
end