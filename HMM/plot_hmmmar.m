%% Load model
load outputs/hmm_33.mat
load outputs/gamma_33.mat

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

% get the MAR coefficients for state k from the estimated model hmm
figure;
tiledlayout(k,order)
for i = 1:k
    W = getMARmodel(hmm,i);
    n_regions = size(W,2);
    start_idx = 1;
    end_idx = n_regions;
    for o = 1:order
        nexttile
        imagesc(W(start_idx:end_idx,:)); colorbar;
        title(sprintf('State %d, AR lag %d',i,o)) % title oaky in this case because timelag==1 (lapse between lags)

        start_idx = end_idx + 1;
        end_idx = end_idx + n_regions;
    end
end
W = getMARmodel(hmm,3);
figure; imagesc(W(1:50,:)); colorbar;
figure; imagesc(W(51:100,:)); colorbar;
figure; imagesc(W(101:150,:)); colorbar;

% get the covariance, correlation and partial matrices for state k, from the estimated model hmm
[covmat,corrmat,icovmat,icorrmat] = getFuncConn(hmm,1);

% also write down modelling choices
% suggestion: cross-validation to try out different options
% also present what people in the paper did
% also consider free energy

A = [4 5; 7 8];
figure; imagesc(A); colorbar;