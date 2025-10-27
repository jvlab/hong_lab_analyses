function status = checkConsist(da,opts)
% JDD - checks for redundancies and consistency across files in a set.
% 10/27/25

numSets = length(da);

for setindx = 1:numSets
    numFiles = length(da{setindx});
    for fileindx = 1:numFiles
        isTarget = da{setindx}{fileindx}.response_amplitude_stim.is_target;
        stim = da{setindx}{fileindx}.response_amplitude_stim.stim(isTarget);
        
        % check these for repeats, and rename these repeats.
        [A,~,ic] = unique(stim);
        if(length(A)<length(stim))
            freqs = accumarray(ic,1)
            idx = find(freqs>1)
            
            for rep = 1:length(idx)
                match_mask = strcmp(A(idx(rep)),stim);
                matchLoc = find(match_mask)
                
            end
            
        end
        
    end
end
    




status=0;
end

