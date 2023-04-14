%% Plot viterbi path
load outputs/vpath_33_full_onpower.mat
figure; 
plot(vpath(1:2000));
title('Viterbi path')
yticks([1 2 3])