% Create a a mixed signal (sinusoid) as a linear combination of individual
% sine waves.
%
% Usage:
%   mixed_signal, sinusoids] = create_sinusoid('fs', 'value', 'T', 'value', 'freqs', 'value')
%
% Inputs:
%   - 'fs'    - [integer] samples per second 
%   - 'T'     - [integer] signal length in seconds
%   - 'freqs' - [cell of integers] frequencies of sinusoids
%
% Output: 
%  - mixed signal - (1,t) array, mixed signal (linear combination of sinusoids)
%  - sinusoids    - (1,n_signals) cell array, individual sinusoids

function [mixed_signal, sinusoids] = create_sinusoid(fs, T, freqs)

    % signal components
    n_signals = size(freqs,2); % number of sinusoids
    dt = 1/fs;                 % seconds per sample
    t = (0:dt:T-dt);           % array of time points
    
    % create sinusoids and add them together
    sinusoids = cell(1,n_signals);
    for i = 1:n_signals
        sinusoids{i} = sin(2*pi*freqs{i}*t);
    end
    mixed_signal = sum(cell2mat(sinusoids'));

end