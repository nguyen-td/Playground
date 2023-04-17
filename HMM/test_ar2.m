%% Generate simulated EEG data from AR model
fs = 200;            % sampling frequency
T = 4;               % trial signal length in seconds
trials = 30;         % number of trials
N = fs * T * trials; % number of samples/data points
M = 62;              % number of channels
P = 10;
K = ceil(M^2/10);

%% Run AR model
[data, Arsig, x, lambdamax] = gen_ar2(M, N, P, K);
outs_ar2 = {data, Arsig, x, lambdamax};
imagesc(cov(x')); colorbar;

histogram(data)
plot(data(1,1:500))
imagesc(Arsig(:,:,5)); colorbar;

%% Get the MAR coefficients for state k from the estimated model hmm
figure;
tiledlayout flow
for i = 1:P
    nexttile
    imagesc(Arsig(:,o)); colorbar;
    title(sprintf('AR lag %d',i)) % title okay in this case because timelag==1 (lapse between lags)
end