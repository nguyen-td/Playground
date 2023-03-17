%% Load model
load outputs/hmm_42.mat
load outputs/gamma_42.mat

%% Plot
% state time courses
k = size(Gamma,2);
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

% get the covariance, correlation and partial matrices for state k, from the estimated model hmm
[covmat,corrmat,icovmat,icorrmat] = getFuncConn(hmm,1);

% get the MAR coefficients for state k from the estimated model hmm
W = getMARmodel(hmm,1);
% use margarita's tile
title('MAR coefficients W')
figure; imagesc(W(1:50,:)); colorbar;
figure; imagesc(W(51:100,:)); colorbar;
figure; imagesc(W(101:150,:)); colorbar;
% also write down modelling choices
% suggestion: cross-validation
% also present what people in the paper did
% also consider free energy

A = [4 5; 7 8];
figure; imagesc(A); colorbar;