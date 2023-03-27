%% Load model
load outputs/hmm_3100_full.mat
% load outputs/hmm_33.mat
load outputs/gamma_3100_full.mat

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

%% Plot state transition matrices
disp(hmm.P)
imagesc(hmm.P); colorbar();
title('State transition matrix')

A = getTransProbs(hmm);
disp(A)
imagesc(A); colorbar;
title('State transition matrix')

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

%% Get noise covariance matrices
% get noise covariance matrix for each state and compare two functions for covtype='full'
figure; 
tiledlayout(2,k)
for mat = 1:2
    if mat == 1
        for i = 1:k
            nexttile
            imagesc(hmm.state(i).Omega.Gam_rate); colorbar;
            title(sprintf('hmm.state(%d).Omega.Gam_rate',i))
        end
    else
        for i = 1:k
            [covmat,~,~,~] = getFuncConn(hmm,i);
            nexttile
            imagesc(covmat); colorbar;
            title(sprintf('covmat of state %d',i))
        end
    end
end

% get noise covariance matrix for each state for covtype='diag'
figure; 
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(diag(hmm.state(i).Omega.Gam_rate)); colorbar;
    title(sprintf('state %d',i))
end

% get shared full noise covariance matrix for covtype='sharedfull'
figure;
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(hmm.Omega.Gam_rate); colorbar;
    title(sprintf('state %d',i))
end

% get noise covariance matrix for each state for covtype='shareddiag'
figure;
tiledlayout(1,k)
for i = 1:k
    nexttile
    imagesc(diag(hmm.Omega.Gam_rate)); colorbar;
    title(sprintf('state %d',i))
end

%% Plot viterbi path
load outputs/vpath_33_diag.mat
figure; 
plot(vpath(1:100));
title('Viterbi path')
yticks([1 2 3])

%% Generate new data
load T
[X_sim,~,~] = simhmmmar(T,hmm);

load X
covtype = 'full'; % adjust
figure;
plot(X(1:500,1))
hold on;
plot(X_sim(1:500,1))
legend('original','simulated')
title(strcat("Original and simulated data, covtype='",covtype,"'"))

%% some notes

% suggestion: cross-validation to try out different options
% also present what people in the paper did
% also consider free energy