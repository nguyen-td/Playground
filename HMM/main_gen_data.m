% Script to generate state-dependent time series data.

%% Add paths
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/plugins/roiconnect/libs/mvgc_v1.0/'))
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Playground/HMM/'))

%% Generate data
% manually define transition matrix
A = [0.98 0.002 0.018;
    0.002 0.997 0.001;
    0.015 0.003 0.9820];

n_states = 3;
P = 10; % number of lags
[data_gt, W, Sigma, sim_gamma] = gen_state_ar(n_states, P, A);

%% Plot transition matrix
figure;
imagesc(A); colorbar;
title('Transition probability matrix A')

%% Save ground truth MVAR model parameters
% DIROUT = 'outputs/'; % change if needed
% if ~exist(DIROUT); mkdir(DIROUT); end
% 
% W_fname = sprintf(strcat(DIROUT, 'W_sim_gt_%d%d.mat'), n_states, P); 
% Sigma_fname = sprintf(strcat(DIROUT, 'Sigma_sim_gt_%d%d.mat'), n_states, P);
% 
% save(W_fname, 'W')
% save(Sigma_fname, 'Sigma')
%     
%% Plots for sanity checks
%     titles = ["MVAR for state 1", "MVAR for state 2", "MVAR for state 3"];
%     plt_all_noisecovs(Sigma, 0, titles)
%     plt_all_Ws(W, P, titles)

%% Save data
% sim_gamma_name = sprintf(strcat(DIROUT, 'Viterbi_sim_gt_%d%d.mat'), n_states, P);
% data_gt_name = sprintf(strcat(DIROUT, 'data_sim_gt_%d%d.mat'), n_states, P);
% 
% save(sim_gamma_name, 'sim_gamma')
% save(data_gt_name, 'data_gt')

%% Plot ground truth time series
% % raw time series
% figure;
% T = (1:size(data_gt, 1)) / fs;
% plot(T, data_gt)
% xlabel('time [s]')
% title(sprintf('Raw ground truth data, %d channels', n_chans))
% 
% % PSDs
% data_gt_pxx = log(abs(pwelch(data_gt, 100, 50)));
% figure;
% tiledlayout(2, P/2)
% for ilag = 1:P
%     nexttile
%     plot(data_gt_pxx(:, ilag))
%     ylabel('log(abs(Power))')
%     title(sprintf('PSDs of generated data, channel %d', ilag))
% end
% 
% % simulated Viterbi path
% figure;
% plot(T, sim_gamma)
% title('Simulated Viterbi path (state time course')
% xlabel('time [s]')
% ylabel('state')
% 
% % %% Create a train_hmmmar compatible file for the 'load_data' parameter
% % out_genar = {data_gt};
% % save(sprintf(strcat(DIROUT,'gen_%d%d_eeg.mat'), n_states, n_chans), 'out_genar')
