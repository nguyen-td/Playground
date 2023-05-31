% Load ground truth data AR data that was previously generated from the gen_ar2.m script called in train_hmmmar.m. 
%
% Input:
%   f_path - file path to .mat file
%
% Outputs:
%   X         - (n_chan, n_obs) AR-generated time series (not multi-trial)
%   Arsig     - (n_chan, n_chan, n_lags) autoregression coefficient matrix
%   x         - (n_chan, n_obs) noise time series
%   lambdamax - largest eigenvalue (see gen_ar2.m)

function [X, Arsig, x, lambdamax] = load_groundtruth_ar(f_path)
    load(f_path)
    X = out_genar{1}';
    Arsig = out_genar{2};
    x = out_genar{3};
    lambdamax = out_genar{4};
end
