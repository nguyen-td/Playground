% Generate state-dependent time series data. Given a transition probability matrix A, an MVAR model generates 
% time points for the given state. Every [t_state] seconds, a new state is determined. A separate ground truth 
% MVAR model is estimated for each state.
%
% Mathematically, the probability that a time point x(t) is generated given the actual state at 
% time t, denoted by m(t), is given by P(x(t) | m(t) = S_j). The probability a_ij of switching 
% from state S_j to state S_i at time t (the actual states at t are given by m(t) and m(t-1) respectively, 
% is given by a_ij = P(m(t) = S_j | m(t-1) = S_i).
%
% Input:
%   n_states - number of hidden states to simulate
%   P        - number of lags
%   A        - (n_states, n_states) transition probability matrix

%
% Outputs:
%   data_gt   - (fs * 4 * 30, n_states) ground truth time series generated by n_states MVAR models 
%   W         - (...) regression coefficient matrices
%   Sigma     - (...) noise covariance matrices
%   sim_gamma - (fs * 4 * 30, n_states) state time course

function [data_gt, W, Sigma, sim_gamma] = gen_state_ar(n_states, P, A)
    % check that the transition matrix is correctly defined
    if ~all(sum(A, 2))
        error('Rows of the transition matrix A must sum up to 1')
    end

    if ~(size(A, 1) && size(A, 2))
        error('The transition matrix A must be a square matrix')
    end

    if ~(size(A, 1) == n_states)
        error('The rows/columns of the transition matrix A must be equal to the number of states')
    end
    
    % simulate AR data using gen_ar2.m
    fs = 200;              % sampling frequency
    t = 4;                 % trial signal length in seconds
    n_trials = 30;         % number of trials
    N = fs * t * n_trials; % number of samples/data points
    M = 10;                % number of channels
    K = ceil(M^2 / 10);
    
    data = cell(1, n_states);
    for istate = 1:n_states
        data{istate} = gen_ar2(M, N, P, K);
    end
    
    % fit MVGC MVAR models to generate "ground truth" MVAR models 
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
    
    % generate individual time series for each state
    data_states = cell(1, n_states);
    for istate = 1:n_states
        X = var_to_tsdata(W{istate}, Sigma{istate}, n_samples, n_trials);
        data_states{istate} = X(:,:)';
    end
    
    % define state probabilities
    states = cell(1, n_states);
    for istate = 1:n_states
        states{istate} = A(istate, :);
    end
    
    % set some parameters
    Pi = [1, 0, 0];           % initial probability distribution, equivalent to options.Pistructure
    t_state = 1;              % duration in seconds of each state segment
    obs_state = fs * t_state; % number of observations/time points during each state segment 
    n_transitions = N / obs_state;
    if ~all(sum(Pi))
        error('The initial probability distribution Pi must sum up to 1')
    end
    
    % generate ground truth time series
    data_gt = zeros(N, n_chans);
    sim_gamma = zeros(1, N); % time course of states (Viterbi path in HMM)
    
    [~, m_t] = max(Pi);
    data_gt(1:obs_state, :) = data_states{m_t}(1:obs_state, :); % first [t_state] seconds
    sim_gamma(1:obs_state) = m_t;
    
    for itrans = 1:n_transitions-1
        for istate = 1:n_states
            if rand < A(m_t, istate)
                m_t = istate; % determine new state
            end
        end
        start_idx = itrans * obs_state + 1;
        end_idx = start_idx + obs_state - 1;
        data_gt(start_idx:end_idx, :) = data_states{m_t}(start_idx:end_idx, :); % sample data points given the current state
        sim_gamma(start_idx:end_idx) = m_t;
    end
end