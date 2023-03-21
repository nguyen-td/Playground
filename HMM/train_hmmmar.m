% Based on run_HMMMAR2.m
% 
% Inputs:
%   k - number of states
%   order - order of the AR model

function train_hmmmar(k, order)
    % addpath
    addpath_hmm

    % prepare data    
    N = 25; % subjects
    Q = 4; % sessions per subject
    ttrial = 500; % time points
    nregions = 50; % regions or voxels
    Y = zeros(N*Q,2); % design matrix with conditions
    X = randn(Q*N*ttrial,nregions); % all data concatenated
    T = ttrial * ones(N*Q,1);  % length of data for each session
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
    
    % set up HMM
    configurations = {}; 
    configurations.order = order;
    configurations.dropstates = 0;
    configurations.verbose = 0;
    configurations.cyc = 50;
    configurations.initcyc = 5;
    configurations.K = k;
    configurations.zeromean = 1; 
    configurations.covtype = 'sharedfull';
            
    % run HMM 
    tic
    [hmm, Gamma] = hmmmar(X,T,configurations); 
    t_end = toc;

    % save model outputs and create logfile to store training time
    DIROUT = 'outputs/'; % change if needed
    if ~exist(DIROUT); mkdir(DIROUT); end

    felapsed_time = fopen(sprintf(strcat(DIROUT,'time_%d%d.txt'), k, order),'w');
    fprintf(felapsed_time,'Elapsed time is %d seconds.',t_end);
    fclose(felapsed_time);

    hmm_name = sprintf(strcat(DIROUT,'hmm_%d%d.mat'), k, order); 
    gamma_name = sprintf(strcat(DIROUT,'gamma_%d%d.mat'), k, order); 
    save(hmm_name, 'hmm', '-v7.3') % saving variables > 2GB
    save(gamma_name, 'Gamma')
end

