function stimulus_list = findRepeatStimuli(S)
% JDD 11/20
% S is a single table.

stimulus_all = S.Properties.RowNames;

stimulus_list = {};


stim_short = stimulus_all;

for stimindx = 1:length(stimulus_all)
    stim_short{stimindx} = stimulus_all{stimindx}(1:end-4);
    if stim_short{stimindx}(1:5)=='diag_'
        stim_short{stimindx} = stim_short{stimindx}(6:end);
    end
end

while length(stim_short)>1
    stim1 = stim_short{1};
    findx = strcmp(stim_short,stim1);
    if(sum(findx)>1)
        stimulus_list = [stimulus_list;stimulus_all(findx)];
    end
    stimulus_all(findx) = [];
    stim_short(findx) = [];    
end
    

    
    

end


