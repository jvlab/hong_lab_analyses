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
        % The dataset indicates which odorants were part of the experiment.
        isTarget = da{setindx}{fileindx}.response_amplitude_stim.is_target;
        stim = da{setindx}{fileindx}.response_amplitude_stim.stim(isTarget);
        
        % check these for repeats, and rename these repeats.
        [A,~,ic] = unique(stim);
        if(length(A)<length(stim))
            status = string('found repeat stimuli');
            % Anything larger than 1 indicated a repeat.    
            freqs = accumarray(ic,1);
            idx = find(freqs>1);
            % This loop for the case where there are multiple, distinct
            % repeats.
            for odor = 1:length(idx)
                match_mask = strcmp(A(idx(odor)),stim);
                matchLoc = find(match_mask);
                % This loop is for the unseen case where an odor is
                % repeated more than once.
                for repeat = 2:length(matchLoc)
                    addpos = find(stim{matchLoc(repeat)} == ' ',1,'first');
                    stim{matchLoc(repeat)} = insertAfter(stim{matchLoc(repeat)},addpos-1,opts.rep_char);                    
                end
            end            
        end
        
        % check the trial info stimuli, if available. 
         
        if(isfield(da{setindx}{fileindx},'trial_info'))
            
            
        end
        
           
        daout{setindx}{fileindx}.response_amplitude_stim.stim = stim;
        
        
    end
end



end

