% Generate state-dependent time series data. According to a transition probability matrix A, 
% an MVAR model generates time points for the given state. Every 5 seconds, a new state is determined. 
%
% Mathematically, the probability b_j that a time point x(t) is generated given the actual state at 
% time t, denoted by m(t), is given by b_j = P(x(t) | m(t) = S_j). The probability a_ij of switching 
% from state S_j to state S_i at time t (the actual states at t are given by m(t) and m(t-1) respectively, 
% is given by a_ij = P(m(t) = S_j | m(t-1) = S_i).
% 
% We define the matrix A and the parameters of the MVAR model manually. Also, we set the number of simulated state to 3. 

%% Add paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/plugins/roiconnect/libs/mvgc_v1.0/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'))
n_states = 3;

%% Define parameters
% transition matrix
A = [0.98 0.002 0.018;
    0.002 0.997 0.001;
    0.015 0.003 0.9820];

if all(sum(A, 2))
    error('Rows of the transition matrix A must sum up to 0')
end

%% Simulate AR data using gen_ar2.m
fs = 200;            % sampling frequency
t = 4;               % trial signal length in seconds
n_trials = 30;         % number of trials
N = fs * t * n_trials; % number of samples/data points
M = 10;              % number of channels
P = 10;              % number of lags
K = ceil(M^2 / 10);

data = cell(1, n_states);
for istate = 1:n_states
    data{istate} = gen_ar2(M, N, P, K);
end

%% Plots for sanity checks
% raw data time series
figure;
t = (1:size(data{1}, 2)) / fs;
for istate = 1:n_states
    hold on;
    plot(t, data{istate})
end
xlabel('time [s]')
legend('time series 1', 'time series 2', 'time series 3')
title('Generated MVAR data (raw data)')

% PSDs
data_pxx = cell(1, n_states);
for istate = 1:n_states
    data_pxx{istate} = log(abs(pwelch(data{istate}', 100, 50)));
end

figure;
tiledlayout(2, P/2)
for ilag = 1:P
    nexttile
    for istate = 1:n_states
        hold on;
        plot(data_pxx{istate}(:, ilag))
    end
    ylabel('log(abs(Power))')
    legend('time series 1', 'time series 2', 'time series 3')
    title(sprintf('PSDs of generated data, channel %d', ilag))
end

%% Fit MVGC MVAR models to generate "ground truth" MVAR models 
n_chans = size(data{1}, 1);
n_samples = size(data{1}, 2) / n_trials;

% reshape time series back to being multi-trial, then fit MVGC MVAR model
data_trial = cell(1, n_states);
W = cell(1, n_states);       % regression coefficient matrices
Sigma = cell(1, n_states);   % noise covariance matrices

for istate = 1:n_states
    data_trial{istate} = reshape(data{istate}, n_chans, n_samples, n_trials); 
    [A, SIG, ~] = tsdata_to_var(data_trial{istate}, P, 'OLS');

    W{istate} = A;
    Sigma{istate} = SIG;
end

%% Save ground truth MVAR model parameters
DIROUT = 'outputs/'; % change if needed
if ~exist(DIROUT); mkdir(DIROUT); end

W_fname = sprintf(strcat(DIROUT, 'W_sim_gt_%d%d.mat'), n_states, P); 
Sigma_fname = sprintf(strcat(DIROUT, 'Sigma_sim_gt_%d%d.mat'), n_states, P);

save(W_fname, 'W')
save(Sigma_fname, 'Sigma')

%% Plots for sanity checks
titles = ["MVAR for state 1", "MVAR for state 2", "MVAR for state 3"];
plt_all_noisecovs(Sigma, 0, titles)
plt_all_Ws(W, P, titles)

%% Generate simulated ground truth state-dependent time series
% first, generate individual time series for each state
data_states = cell(1, n_states);
for istate = 1:n_states
    X = var_to_tsdata(W{istate}, Sigma{istate}, n_samples, n_trials);
    data_states{istate} = X(:,:)';
end