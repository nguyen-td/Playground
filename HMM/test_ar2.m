% generate simulated EEG data from AR model
M = 62;
N = 10000;
P = 10;
K = ceil(M^2/10);

[data, Arsig, x, lambdamax] = gen_ar2(M, N, P, K);
outs_ar2 = {data, Arsig, x, lambdamax};
imagesc(cov(x')); colorbar;
histogram(data)