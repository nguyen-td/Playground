%% Load model
load outputs/hmm_33_full_onpower.mat
% load outputs/hmm_33.mat
load outputs/gamma_33_full_onpower.mat

%% Get parameters
% number of states and MAR order
k = size(Gamma,2);
order = size(hmm.state,2);

%% Get the MAR coefficients for state k from the estimated model hmm
figure;
tiledlayout(k,order)
for i = 1:k
    W = getMARmodel(hmm,i);
%     W = hmm.state(i).W.Mu_W; % equivalent
    n_regions = size(W,2);
    start_idx = 1;
    end_idx = n_regions;
    for o = 1:order
        nexttile
        imagesc(W(start_idx:end_idx,:)); colorbar;
        title(sprintf('State %d, AR lag %d',i,o)) % title okay in this case because timelag==1 (lapse between lags)

        start_idx = end_idx + 1;
        end_idx = end_idx + n_regions;
    end
end