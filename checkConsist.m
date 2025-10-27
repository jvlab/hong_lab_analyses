function [status,daout,optsout] = checkConsist(da,opts)
% JDD - checks for redundancies and consistency across files in a set.
% 10/27/25
if(nargin<=1)
    opts = struct;
end

if(~isfield(opts,'rep_char'))
    opts.rep_char = '+';
end

optsout = opts;
numSets = length(da);

daout = da;

status = string('consistent within sets, no repeats found');

for setindx = 1:numSets
    numFiles = length(da{setindx});
    for fileindx = 1:numFiles
        isTarget = da{setindx}{fileindx}.response_amplitude_stim.is_target;
        stim = da{setindx}{fileindx}.response_amplitude_stim.stim(isTarget);
        
        % check these for repeats, and rename these repeats.
        [A,~,ic] = unique(stim);
        if(length(A)<length(stim))
            status = string('found repeat stimuli');
            freqs = accumarray(ic,1);
            idx = find(freqs>1);
            
            for rep = 1:length(idx)
                match_mask = strcmp(A(idx(rep)),stim);
                matchLoc = find(match_mask);
                
                for rep = 2:length(matchLoc)
                    % add the character at the first space
                    
                    addpos = find(stim{matchLoc(rep)} == ' ',1,'first');
                    stim{matchLoc(rep)} = insertAfter(stim{matchLoc(rep)},addpos-1,opts.rep_char);
                    
                end
            end            
        end
        daout{setindx}{fileindx}.response_amplitude_stim.stim = stim;
        
    end
end



end

