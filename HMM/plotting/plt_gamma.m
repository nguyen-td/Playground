%% Load model
load ../outputs/hmm_33_full_onpower.mat
load ../outputs/gamma_33_full_onpower.mat

%% Get parameters
% number of states and MAR order
k = size(Gamma,2);
order = size(hmm.state,2);

%% Plot state time courses (Gamma)
figure; imagesc(Gamma'); colorbar();
title(sprintf('State time course (%d states)', k))

figure;
tiledlayout(k,1)
for i = 1:k
    nexttile
    plot(Gamma(1:100,i));
    title(sprintf('State %d',i))
end