function [status,tableList,daout,optsout] = checkConsist(da,opts)
% JDD - checks for redundancies and consistency across files in a set.
% 10/27/25
% Duplicate stimuli are renamed according to opt.rep_char and treated as
% distinct.
%
% Checks that the glomeruli are the same across data sets.
%
% Checks that the stimuli are the same for files within a set.
%
% Checks that the trial stimuli are consistent.
%
% Outputs a data table with glomeruli columns and stimuli rows.
% Outputs a data structure, very similar to the file structure. 
if(nargin<=1)
    opts = struct;
end

if(~isfield(opts,'kdf'))
    opts.kdf = {'max_peak','mean_peak'};
end
if(~isfield(opts,'rep_char'))
    opts.rep_char = '+';
end

optsout = opts;
numSets = length(da);




daout = da;

status = string('beginning ');

% Check that the trial stimuli are consistent with the averaged stimuli

tableList = cell(numSets,1);

% Need to check that the glomeruli are consistent.

glomeruli_base = da{1}{1}.rois.glomeruli;

for setindx = 1:numSets
    numFiles = length(da{setindx});
    tableList{setindx} = cell(numFiles,1);
    for fileindx = 1:numFiles
        % Check to get the field name of the data.
        for field_name = opts.kdf
            if(isfield(da{setindx}{fileindx}.response_amplitude_stim,field_name{1}))
                datafield = field_name{1};
            end
        end
                
        % The dataset indicates which odorants were part of the experiment.
        isTarget = da{setindx}{fileindx}.response_amplitude_stim.is_target;
        stim = da{setindx}{fileindx}.response_amplitude_stim.stim(isTarget);
        trialTarget = da{setindx}{fileindx}.trial_info.is_target;
        stim_trial = da{setindx}{fileindx}.trial_info.stim(trialTarget);
        
        for stim_val = 1:length(stim)
            for trial_val = 1:opts.trial_repeats
                if(~strcmp(stim(stim_val),stim_trial((stim_val-1)*opts.trial_repeats+trial_val)))    
                   error('response stimuli and trial stimuli do not line up');
                end
            end
        end
        % check these for repeats, and rename these repeats.
        [A,~,ic] = unique(stim);
        if(length(A)<length(stim))
            status = status + string('found repeat stimuli in');
            % Anything larger than 1 indicated a repeat.    
            freqs = accumarray(ic,1);
            idx = find(freqs>1);
            % This loop for the case where there are multiple, distinct
            % repeats.
            for odor = 1:length(idx)
                match_mask = strcmp(A(idx(odor)),stim);
                matchLoc = find(match_mask);
                % This loop is for the yet unseen case where an odor is
                % repeated more than once.
                for repeat = 2:length(matchLoc)
                    addpos = find(stim{matchLoc(repeat)} == ' ',1,'first');
                    stim{matchLoc(repeat)} = insertAfter(stim{matchLoc(repeat)},addpos-1,repmat(opts.rep_char,1,repeat-1));                    
                    for trial=1:opts.trial_repeats
                        stim_trial{(matchLoc(repeat)-1)*opts.trial_repeats+trial} = stim{matchLoc(repeat)};
                    end
                
                end
            end            
        end
        
        % check the trial info stimuli, if available. 
        
        % I know it is, because the script didn't exit earlier.
        if(isfield(da{setindx}{fileindx},'trial_info'))            
        end
        
        % bookkeeping.
        daout{setindx}{fileindx}.response_amplitude_stim.stim = stim;
        daout{setindx}{fileindx}.trial_info.stim = stim_trial;
        daout{setindx}{fileindx}.response_amplitude_stim.panel(~isTarget) = [];
        daout{setindx}{fileindx}.trial_info.panel(~trialTarget) = [];

        % I haven't remove anything, nor have I added anything or moved anything. 
        daout{setindx}{fileindx}.response_amplitude_stim.(datafield)(~isTarget,:)=[];
        daout{setindx}{fileindx}.response_amplitude_trials.(datafield)(~trialTarget,:)=[];
        T = array2table(daout{setindx}{fileindx}.response_amplitude_stim.(datafield));
        glomeruli = daout{setindx}{fileindx}.rois.glomeruli;
        if(~isequal(glomeruli_base,glomeruli))
            errmsg = string('Mismatched glomeruli is set ') + string(setindx) + string(' and file ') + string(fileindx);
            error(errmsg)
        end
        for indx=1:length(glomeruli)
            glomeruli{indx} = strrep(glomeruli{indx},'/','_');
        end
    
        T.Properties.VariableNames = glomeruli;
        T.Properties.RowNames = stim;
        
        tableList{setindx}{fileindx} = T;
    end
    
    % Compare the files' stimulus sets with one another.
    % We will use the first file as the reference stim set
    for fid = 2:numFiles
        set1 = daout{setindx}{1}.response_amplitude_stim.stim;
        set2 = daout{setindx}{fid}.response_amplitude_stim.stim;
        if(length(unique([set1,set2])) > min(length(set1),length(set2)))
            error('Files within a set must have the same stimuli');
        end
    end
    % If we made it out of the loop we have a consistent stimulus set. 
    stimlist = daout{setindx}{1}.response_amplitude_stim.stim;
    
    % Cycle through the stimuli
    for st = 1:length(stimlist)
       % Identify the number of components
       if(sum(stimlist{st}=='@') > 1)
            splitstim = split(stimlist{st});
            odorString = string();
            concString = string();
            for rep=1:4:length(splitstim)
                odorString = odorString + '+' + splitstim(rep);
                concString = concString + ',' + string(splitstim(rep+2));
            end
            odorString = strip(odorString,'+');
            concString = strip(concString,',');
            concString = string('(')+concString+string(')');
       
       else
           splitstim = split(stimlist{st});
           odorString = splitstim(1);
           concString = string('(') + splitstim(3) + string(')');
           
       end
       
       stimlist{st} = char(odorString + concString + string(opts.set_names(setindx)));
        % incorporate the new stimulus names in the data structure and the
        % tables.
       for fileindx = 1:numFiles
           daout{setindx}{fileindx}.response_amplitude_stim.stim = stimlist;
           T = tableList{setindx}{fileindx};
           T.Properties.RowNames = stimlist; % This seems to be the only way this works. 
           tableList{setindx}{fileindx} = T;
       end
       
    end
    

end

status = status + string('Success'); 

end

