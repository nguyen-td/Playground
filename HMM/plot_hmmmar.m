%% Load model
load outputs/hmm_33_diag.mat
% load outputs/hmm_33.mat
load outputs/gamma_33_sharedfull.mat

%% Plot
% number of states and MAR order
k = size(Gamma,2);
order = size(hmm.state,2);

% state time courses
figure; imagesc(Gamma'); colorbar();
title(sprintf('State time course (%d states)', k))

figure;
tiledlayout(k,1)
for i = 1:k
    nexttile
    plot(Gamma(:,i));
    title(sprintf('State %d',i))
end

% state transition matrix
disp(hmm.P)
imagesc(hmm.P); colorbar();
title('State transition matrix')

A = getTransProbs(hmm);
disp(A)
imagesc(A); colorbar;
title('State transition matrix')

% get the MAR coefficients for state k from the estimated model hmm
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

% get noise covariance matrix for each state and compare two functions
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

% get noise covariance matrix for each state
figure; imagesc(hmm.state(1).Omega.Gam_rate); colorbar;

% get shared full noise covariance matrix
Omega_33_sharedfull = hmm.Omega;
save('Omega_33_sharedfull','Omega_33_sharedfull')

% get viterbi path
load X
load T
tic
[viterbipath] = hmmdecode(X,T,hmm,1); % 1350 seconds
toc
save('viterbipath','viterbipath')
plot(viterbipath)


% suggestion: cross-validation to try out different options
% also present what people in the paper did
% also consider free energy