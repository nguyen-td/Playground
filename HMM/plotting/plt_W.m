%% Load model
load outputs/hmm_110_diag_eeg.mat
% load outputs/hmm_33.mat
load outputs/gamma_110_diag_eeg.mat

%% Get parameters
% number of states and MAR order
k = size(Gamma,2);
order = size(hmm.state(1).Omega.Gam_irate, 1);

%% Get the MAR coefficients for state k from the estimated model hmm
figure;
% tiledlayout(k,order)
tiledlayout(2, 5)
for i = 1:k
    W = getMARmodel(hmm,i);
%     W = hmm.state(i).W.Mu_W; % equivalent
    n_regions = size(W,2);
    start_idx = 1;
    end_idx = n_regions;
    % min_val = min(W, [], 'all');
    % max_val = max(W, [], 'all');
    for o = 1:order
        nexttile
        imagesc(W(start_idx:end_idx,:)); colorbar; clim([min_val max_val]);
        title(sprintf('State %d, AR lag %d',i,o)) % title okay in this case because timelag==1 (lapse between lags)

        start_idx = end_idx + 1;
        end_idx = end_idx + n_regions;
    end
end

%% Set up AR model and generate ground truth data
fs = 200;            % sampling frequency
t = 4;               % trial signal length in seconds
trials = 30;         % number of trials
N = fs * t * trials; % number of samples/data points
M = 10;              % number of channels
P = 10;              % number of lags
K = ceil(M^2 / 10);

[data, Arsig, x, lambdamax] = gen_ar2(M, N, P, K);

%% Get the MAR coefficients the (original) AR model
figure;
tiledlayout(2, 5) % to be adjusted
for i = 1:P
    nexttile
    imagesc(Arsig(:,:,i)); colorbar; clim([min_val max_val]);
    title(sprintf('AR lag %d',i)) % title okay in this case because timelag==1 (lapse between lags)
end

%% Compute Eudlidian norm
W_empty = ones(P, P);
tiledlayout(5, 2) % to be adjusted
start_idx = 1;
end_idx = n_regions;
for i = 1:P
    V = squeeze(Arsig(:,:,i) - W(start_idx:end_idx,:));
    l2_norm = sqrt(V(:)' * V(:));

    nexttile
    imagesc(W_empty); colormap(cm18);
    text(5, 5, sprintf('d = %d', l2_norm), 'VerticalAlignment', 'middle', 'VerticalAlignment', 'middle');
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    title(sprintf('L2 norm, AR lag %d',i)) 

    start_idx = end_idx + 1;
    end_idx = end_idx + n_regions;
end