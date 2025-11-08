% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

% dependencies : fileToRaw.m, checkConsist.m, lookForSetWideHoles,
% fillInNaNs, makePlots_1, calcResp, makePlots_quantile.

opts = struct;
opts.rep_char = '+';
opts.trial_repeats = 3;
opts.set_names = {'kiTC','meTC','vaTC'};
opts.kdf = {'max_peak','mean_peak'}; % Known data fields.
opts.suppressoutput = true;
opts.interactive = false;
opts.restore_size = true;
opts.submean = false;
opts.hist_quantiles = [0.05 .25 .5 .75 .96];
opts.hist_bins = 50;

hlid_setup;
% I am assuming that the files that are part of a set are in a separate
% directory. There is not name checking, if a data file is in the folder
% the code will attempt to load it into the set.
Sraw{1} = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Sraw{2} = fileToRaw('../orn_terminals_Oct25/megamat17');
%Sraw{3} = fileToRaw('../orn_terminals_Oct25/monat');
Sraw{3} = fileToRaw('../orn_terminals_Oct25/validation2');

% Check stimulus consistency across files within a set,
% and glomerus consistency across sets.
[~,Sall] = checkConsist(Sraw,opts);

% First, I want to check that each stimulus has a non-NaN value somewhere
% within a set. If an all NaN is found, that stimulus is removed.
% The presence of the glomeruli across sets is determined. the glomeruli
% used appear in a minimum number of files. 
Strimmed = lookForSetWideHoles(Sall,opts);

% Fill in the remaining holes.
% The calls the afalwt interpolator.
[Sfilled,afalwt_fit] = fillInNaNs(Strimmed,opts);

% Generate the first set of plots (raw - trimmed - filled)
%makePlots_1(Sall,Strimmed,Sfilled);

%
% create the resps_set table
[resps_set,resp_range] = calcResp(Sfilled,afalwt_fit,opts);

% Generate the quantile plots
%makePlots_quantile(resps_set,resp_range,opts);

% Merge the sets into an intersection and union of glomeruli.
merged_data = mergeSets(resps_set);

% All good to here. The merged sets are in the merged data cell array.

numSets = length(Sfilled);
dsid = cell(numSets,1);
for setindx = 1:numSets
    numFiles = length(Sfilled{setindx});
    for fileindx = 1:numFiles
        dsid_this = Sraw{setindx}{fileindx}.meta.title;
        dsid_this=strrep(dsid_this,'/','-');
        dsid_this=strrep(dsid_this,'\','-');
        dsid_this=strrep(dsid_this,'_','-');
        while contains(dsid_this,'--')
            dsid_this=strrep(dsid_this,'--','-');
        end
        %
        dsid{setindx}=[dsid{setindx};dsid_this];
    end
end
f_base=struct;
for setindx=1:numSets
    f_base.metadata{setindx}=Sraw{setindx}{1}.meta; %original metadata from Hong Lab
    f_base.dsid{setindx}=dsid{setindx}; %data set ID, with special chars turned into -
end
f_base.resps=resps_set; %original responses
if ~opts.restore_size
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in'; %original field for responses from Hong Lab
else
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in and restoring to original size';
end

% This is just to get by. I need to figure out why there are four different
% lists of stimuli.
stimulusNames = merged_data{1}.Properties.RowNames;
f_base.stimulus_names_orig=char(stimulusNames); %original stimulus names
f_base.stim_labels_orig=char(stimulusNames); %shortened names for plotting
f_base.stimulus_names=f_base.stimulus_names_orig;
f_base.stim_labels=f_base.stim_labels_orig;

