% Generate state-dependent time series data. Given a transition probability matrix A, an MVAR model generates 
% time points for the given state. Every [t_state] seconds, a new state is determined. A separate ground truth 
% MVAR model is estimated for each state.
%
% Mathematically, the probability that a time point x(t) is generated given the actual state at 
% time t, denoted by m(t), is given by P(x(t) | m(t) = S_j). The probability a_ij of switching 
% from state S_j to state S_i at time t (the actual states at t are given by m(t) and m(t-1) respectively, 
% is given by a_ij = P(m(t) = S_j | m(t-1) = S_i).
% 
% We define the matrix A manually. Also, we set the number of simulated state to 3. 

%% Add paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/plugins/roiconnect/libs/mvgc_v1.0/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'))
n_states = 3;

%% Define parameters
% transition matrix
A = [0.98 0.002 0.018;
    0.002 0.997 0.001;
    0.015 0.003 0.9820];

if ~all(sum(A, 2))
    error('Rows of the transition matrix A must sum up to 1')
end

% plot transition probability matrix
figure;
imagesc(A)
title('Transition probability matrix A')

%% Simulate AR data using gen_ar2.m
fs = 200;              % sampling frequency
t = 4;                 % trial signal length in seconds
n_trials = 30;         % number of trials
N = fs * t * n_trials; % number of samples/data points
M = 10;                % number of channels
P = 10;                % number of lags
K = ceil(M^2 / 10);

data = cell(1, n_states);
for istate = 1:n_states
    data{istate} = gen_ar2(M, N, P, K);
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
    [coeff_mat, SIG, ~] = tsdata_to_var(data_trial{istate}, P, 'OLS');

    W{istate} = coeff_mat;
    Sigma{istate} = SIG;
end

%% Save ground truth MVAR model parameters
% DIROUT = 'outputs/'; % change if needed
% if ~exist(DIROUT); mkdir(DIROUT); end
% 
% W_fname = sprintf(strcat(DIROUT, 'W_sim_gt_%d%d.mat'), n_states, P); 
% Sigma_fname = sprintf(strcat(DIROUT, 'Sigma_sim_gt_%d%d.mat'), n_states, P);
% 
% save(W_fname, 'W')
% save(Sigma_fname, 'Sigma')

%% Plots for sanity checks
titles = ["MVAR for state 1", "MVAR for state 2", "MVAR for state 3"];
plt_all_noisecovs(Sigma, 0, titles)
plt_all_Ws(W, P, titles)

%% Generate simulated ground truth state-dependent time series
% generate individual time series for each state
data_states = cell(1, n_states);
for istate = 1:n_states
    X = var_to_tsdata(W{istate}, Sigma{istate}, n_samples, n_trials);
    data_states{istate} = X(:,:)';
end

% define states and store the probabilities
states = cell(1, n_states);
for istate = 1:n_states
    states{istate} = A(istate, :);
end

% set some parameters
Pi = [1, 0, 0];           % initial probability distribution, equivalent to options.Pistructure?
t_state = 1;              % duration in seconds of each state segment
obs_state = fs * t_state; % number of observations/time points during each state segment 
n_transitions = N / obs_state;
if ~all(sum(Pi))
    error('The initial probability distribution Pi must sum up to 1')
end

% generate ground truth time series
data_gt = zeros(N, n_chans);
sim_viterbi = zeros(1, N); % time course of states (Viterbi path in HMM)

[~, m_t] = max(Pi);
data_gt(1:obs_state, :) = data_states{m_t}(1:obs_state, :); % first [t_state] seconds
sim_viterbi(1:obs_state) = m_t;

for itrans = 1:n_transitions-1
    for istate = 1:n_states
        if rand < A(m_t, istate)
            m_t = istate; % determine new state
        end
    end
    start_idx = itrans * obs_state + 1;
    end_idx = start_idx + obs_state-1;
    data_gt(start_idx:end_idx, :) = data_states{m_t}(start_idx:end_idx, :); % sample data points given the current state
    sim_viterbi(start_idx:end_idx) = m_t;
end

% % save data
% sim_viterbi_name = sprintf(strcat(DIROUT, 'Viterbi_sim_gt_%d%d.mat'), n_states, P);
% data_gt_name = sprintf(strcat(DIROUT, 'data_sim_gt_%d%d.mat'), n_states, P);
% 
% save(sim_viterbi_name, 'sim_viterbi')
% save(data_gt_name, 'data_gt')
%% Plot ground truth time series
% raw time series
figure;
T = (1:size(data_gt, 1)) / fs;
plot(T, data_gt)
xlabel('time [s]')
title(sprintf('Raw ground truth data, %d channels', n_chans))

% PSDs
data_gt_pxx = log(abs(pwelch(data_gt, 100, 50)));
figure;
tiledlayout(2, P/2)
for ilag = 1:P
    nexttile
    plot(data_gt_pxx(:, ilag))
    ylabel('log(abs(Power))')
    title(sprintf('PSDs of generated data, channel %d', ilag))
end

% simulated Viterbi path
figure;
plot(T, sim_viterbi)
title('Simulated Viterbi path (state time course')
xlabel('time [s]')
ylabel('state')
