%% Load model
load ../outputs/hmm_33_full_onpower.mat
load ../outputs/gamma_33_full_onpower.mat

%% Get parameters
% number of states and MAR order
k = size(Gamma,2);
order = size(hmm.state,2);

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