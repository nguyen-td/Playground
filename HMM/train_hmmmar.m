% Based on run_HMMMAR_2.m
% 
% Inputs:
%   k        - number of states
%   order    - order of the AR model
%   viterbi  - specify if viterbi path should be saved, boolean (1 or 0)
%   data_mod - data modality, string ('eeg' or 'fmri')

function train_hmmmar(k, order, covtype, viterbi, data_mod)
    % addpath
    addpath_hmm
    
    % generate simulated data
    if strcmpi(data_mod,'fmri')
        [X, T] = gen_fmri;
        out_genfmri = {X, T};
    else
        M = 62; N = 10000; P = 10; K = ceil(M^2/10); % AR parameters
        [data, Arsig, x, lambdamax] = gen_ar2(M, N, P, K);
        out_genar = {data, Arsig, x, lambdamax};
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
    % configurations.exptimelag = 2;
    % configurations.onpower = 1; % Hilbert transform as done in Baker et al (2014), eLife.
    configurations.covtype = covtype;
    % choose P = 100
            
    % run HMM 
    if viterbi
        tic
        [hmm, Gamma, ~, vpath] = hmmmar(X, T, configurations); 
        t_end = toc;
    else
        tic
        [hmm, Gamma] = hmmmar(X, T, configurations); 
        t_end = toc;
    end

    % save model outputs and create logfile to store training time
    DIROUT = 'outputs/'; % change if needed
    if ~exist(DIROUT); mkdir(DIROUT); end

    felapsed_time = fopen(sprintf(strcat(DIROUT,'time_%d%d_',covtype,'.txt'), k, order),'w');
    fprintf(felapsed_time,'Elapsed time is %d seconds.',t_end);
    fclose(felapsed_time);
    
    hmm_name = sprintf(strcat(DIROUT,'hmm_%d%d_',covtype,'.mat'), k, order); 
    gamma_name = sprintf(strcat(DIROUT,'gamma_%d%d_',covtype,'.mat'), k, order); 
    save(hmm_name, 'hmm', '-v7.3') % saving variables > 2GB
    save(gamma_name, 'Gamma')
    if viterbi
        vpath_name = sprintf(strcat(DIROUT,'vpath_%d%d_',covtype,'.mat'), k, order); 
        save(vpath_name, 'vpath') 
    end
    if strcmpi(data_mod,'fmri')
        genfmri_name = strcat(DIROUT,'gen_fmri.mat');
        save(genfmri_name, 'out_genfmri')
    else
        genar_name = strcat(DIROUT,'gen_ar');
        save(genar_name, 'out_genar')
    end
end

