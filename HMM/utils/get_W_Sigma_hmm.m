% Get regression coefficients and covariance matrices from HMM-MAR trained with both parameter values:
% covtype = 'full' and covtype = 'diag'. Outputs are in this order.
% 
% Inputs:
%   k     - number of states
%   order - order of the AR model
%
% Outputs:
%   Ws      - (1, 2) cell array containing (n_chan * n_chan, n_lags) regression coefficient matrices
%   covmats - (1, 2) cell array containing (n_chan, n_chan) noise covariance matrices

function [Ws, covmats] = get_W_Sigma_hmm(k, order)
    covtypes = {'full', 'diag'};
    Ws = cell(1, 2);
    covmats = cell(1, 2);
    
    for icov = 1:length(covtypes)
        load(sprintf(strcat('outputs/hmm_%d%d_', covtypes{icov}, '_eeg.mat'), k, order))
        load(sprintf(strcat('outputs/Gamma_%d%d_', covtypes{icov}, '_eeg.mat'), k, order))
        W = getMARmodel(hmm, 1); % get regression coefficient
        Ws{icov} = W;
        
        if strcmpi(covtypes{icov}, 'full')
            [covmat, ~, ~, ~] = getFuncConn(hmm, k); % noise covariance matrix
            covmats{icov} = covmat;
        else % strcmpi(icov, 'diag')
            covmat = diag(hmm.state(k).Omega.Gam_rate) / (hmm.state(k).Omega.Gam_shape-1); % noise covariance matrix
            covmats{icov} = covmat;
        end
    end
end
