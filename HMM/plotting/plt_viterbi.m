%% Plot viterbi path
load outputs/vpath_310_full_eeg.mat
figure; 
plot(vpath);
title('Viterbi path')
yticks([1 2 3])

% plot histogram