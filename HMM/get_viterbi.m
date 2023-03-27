% Inputs:
%   f_hmm - file name of the HMM model, e.g. hmm_41

function get_viterbi(f_hmm)
    DIROUT = 'outputs/'; % change if needed
    load X
    load T
    load(strcat(DIROUT,f_hmm))

    [viterbipath] = hmmdecode(X,T,hmm,1); % get viterbi path
    save(strcat(DIROUT,'viterbipath'),'viterbipath')
end