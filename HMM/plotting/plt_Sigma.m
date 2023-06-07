%% Load model
cd '/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM'
load outputs/hmm_310_full_eeg.mat
load outputs/gamma_310_full_eeg.mat

%% Get parameters
% number of states and MAR order
k = size(Gamma, 2);
order = size(hmm.state(1).Omega.Gam_irate, 1);

%% Set up AR model and generate ground truth data
fs = 200;            % sampling frequency
t = 4;               % trial signal length in seconds
trials = 30;         % number of trials
N = fs * t * trials; % number of samples/data points
M = 10;              % number of channels
P = 10;              % number of lags
K = ceil(M^2 / 10);

[data, Arsig, x, lambdamax] = gen_ar2(M, N, P, K);

%% Plot noise covariance matric of (original) AR model
covx = cov(x');
imagesc(covx); colorbar;

%% Get noise covariance matrices of estimated HMM-MAR
% get noise covariance matrix for each state for covtype='full'
[covmat, ~, ~, ~] = getFuncConn(hmm,k);
figure; 
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(covmat); colorbar;
    title(sprintf('state %d',i))
end

% get noise covariance matrix for each state for covtype='diag'
figure; 
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(diag(hmm.state(i).Omega.Gam_rate) / (hmm.state(k).Omega.Gam_shape-1)); colorbar;
    title(sprintf('state %d',i))
end

% get shared full noise covariance matrix for covtype='sharedfull'
figure;
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(hmm.Omega.Gam_rate); colorbar;
    title(sprintf('state %d',i))
end

% get noise covariance matrix for each state for covtype='shareddiag'
figure;
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(diag(hmm.Omega.Gam_rate)); colorbar;
    title(sprintf('state %d',i))
end