%% Load model
load outputs/hmm_33.mat
load outputs/gamma_33.mat

%% Plot
% state time courses
imagesc(Gamma'); colorbar();
title('State time course (4 states)')

Gamma_trans = Gamma';
plot(Gamma_trans(1,1:100)');

% state transition matrix
state_trans_mat = hmm.P;
disp(state_trans_mat)
imagesc(state_trans_mat); colorbar();
title('State transition matrix')

% get the covariance, correlation and partial matrices for state k, from the estimated model hmm
[covmat,corrmat,icovmat,icorrmat] = getFuncConn(hmm,1);

% get the MAR coefficients for state k from the estimated model hmm
W = getMARmodel(hmm,1);

imagesc(covmat); colorbar();

a = 4;
b = 5;
tic
c = a + b;
t_end = toc