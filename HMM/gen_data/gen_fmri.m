% Generate simulated fMRI data
function [X, T] = gen_fmri
    N = 25;                           % subjects
    Q = 4;                            % sessions per subject
    ttrial = 500;                     % time points
    nregions = 50;                    % regions or voxels
    Y = zeros(N*Q,2);                 % design matrix with conditions
    X = randn(Q*N*ttrial,nregions);   % all data concatenated
    T = ttrial * ones(N*Q,1);         % length of data for each session
    Tsubject = Q*ttrial * ones(N,1);  % length of data for each subject
    
    % functional connectivity profile for condition 1
    FC_condition1 = randn(nregions); 
    FC_condition1 = FC_condition1' * FC_condition1;
    % functional connectivity profile for condition 2
    FC_condition2 = randn(nregions); 
    FC_condition2 = FC_condition2' * FC_condition2;
    % note that there is no difference in the mean activity between conditions
    
    for j1 = 1:N
        for j2 = 1:Q
            t = (1:ttrial) + (j2-1)*ttrial + Q*(j1-1)*ttrial;
            n = (j1-1)*Q + j2;
            if rem(j2,2)
                Y(n,1) = 1;
                X(t,:) = X(t,:) * FC_condition1;
            else
                Y(n,2) = 1;
                X(t,:) = X(t,:) * FC_condition2;
            end
            for i = 1:nregions 
                X(t,i) = smooth(X(t,i));
            end
        end
    end
end