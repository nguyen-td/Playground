%% Add paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'))

%% Generate simulated EEG data from AR model
fs = 200;            % sampling frequency
T = 4;               % trial signal length in seconds
trials = 30;         % number of trials
N = fs * T * trials; % number of samples/data points
M = 10;              % number of channels
P = 3;
K = ceil(M^2/10);

%% Run AR model
[data, Arsig, x, lambdamax] = gen_ar2(M, N, P, K);
covx = cov(x');
imagesc(covx); colorbar;
outs_ar2 = {data, Arsig, x, lambdamax};

histogram(x)

%% Get the MAR coefficients from the (original) AR model hmm
figure;
tiledlayout(floor(P/2),2)
for i = 1:P
    nexttile
    imagesc(Arsig(:,:,i)); colorbar;
    title(sprintf('AR lag %d',i)) % title okay in this case because timelag==1 (lapse between lags)
end