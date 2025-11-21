function stimulus_list = findRepeatStimuli(S)
% JDD 11/20

numSets = length(S);

stimulus_all = {};

for setindx = 1:numSets
    stimulus_all = [stimulus_all; S{setindx}{1}.Properties.RowNames];
end
ss = {};

% Strip stimuli of ids
stim_short = stimulus_all;

for stimindx = 1:length(stimulus_all)
    stim_short{stimindx} = stimulus_all{stimindx}(1:end-4);
    if stim_short{stimindx}(1:4)=='diag'
        stim_short{stimindx} = stim_short{stimindx}(5:end);
    end
end

for stimindx = 1:length(stimulus_all)
    stim1 = stim_short{stimindx};
    if strcmp(stim1,string('n'))
        continue
    end
    
    
    ff=strcmp(stim_short,stim1);
    if(sum(ff)>1)
        ss=[ss;stimulus_all(find(ff))]
        stim_short(find(ff))=string('n');
    end
    
    

end


