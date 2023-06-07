%% Load model
load outputs/hmm_310_diag_eeg.mat
% load outputs/hmm_33.mat
load outputs/gamma_310_diag_eeg.mat

%% Plot state transition matrices
disp(hmm.P)
figure;
imagesc(hmm.P); colorbar();
title('State transition matrix')

figure,
A = getTransProbs(hmm);
disp(A)
imagesc(A); colorbar;
title('State transition matrix')