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
makePlots_1(Sall,Strimmed,Sfilled);

%
% create the resps_set table
[resps_set,resp_range] = calcResp(Sfilled,afalwt_fit,opts);

% Generate the quantile plots
makePlots_quantile(resps_set,resp_range,opts);

% Since the odorants within a set are unique and each set is uniquely
% tagged, we don't need to do all of the uniqueness checks on the stimuli.
% I do need to assemble the odorants into a single usnique structure.
numSets = length(Sfilled);
stimulus_list = {};
for setindx = 1:numSets
    stimulus_list = [stimulus_list; Sfilled{setindx}{1}.Properties.RowNames];
end

f_base = struct;
for setindx = 1:numSets
    f_base.metadata{setindx} = Sraw{setindx}{1}.meta;
end

f_base.resps=resps_set; %original responses
if (~opts.restore_size)
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in'; %original field for responses from Hong Lab
else
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in and restoring to original size';
end

ncombm=2; %number of combination modes: 1: intersection, 2: union, set missing to 0
resps_combined=cell(1,ncombm); %combined responses, for saving
resps_concat=cell(1,ncombm); %concatenated responses, for plotting
f_combm=cell(1,ncombm);
s_diag_all_combm=cell(1,ncombm);
u_full_combm=cell(1,ncombm);
v_full_combm=cell(1,ncombm);
s_full_combm=cell(1,ncombm);
coords_all_combm=cell(1,ncombm);
%
s_diag_set=cell(numSets,ncombm);
u_full_set=cell(numSets,ncombm);
v_full_set=cell(numSets,ncombm);
s_full_set=cell(numSets,ncombm);

nglomeruli = length(Sall{1}{1}.Properties.VariableNames);

glomeruli_combined=cell(1,ncombm);
glomeruli_combined{1}=Sall{1}{1}.Properties.VariableNames;
glomeruli_combined{2}=[];
for setindx=1:numSets
    glom_use_set = Sfilled{setindx}{1}.Properties.VariableNames;
    glomeruli_combined{1}=intersect(glomeruli_combined{1},glom_use_set);
    glomeruli_combined{2}=union(glomeruli_combined{2},glom_use_set);
end


if_ok=0;
while (if_ok==0)
    if(opts.interactive)
        if_restrict_stims=getinp('1 to restrict stimuli to a subset before computing pcs of combined sets','d',[0 1]);
    else
        if_restrict_stims=0;
    end
        
    if if_restrict_stims
        disp('**********'); %code from psg_coord_pipe_proc
        disp('Enter a file for the purposes of defining a stimulus set.')
        [sets_tn,ds_tn,sas_tn]=psg_get_coordsets(opts_read,[],[],1);
        subset_typenames=sas_tn{1}.typenames;
        disp('Stimulus set retrieved.');
        disp(subset_typenames');
        disp('**********');
        stims_combined_keep=[];
        for istim=1:length(stimulus_list)
            if length(strmatch(stimiulus_list{istim},subset_typenames,'exact'))==1
                stims_combined_keep(end+1)=istim;
            end
        end   
        disp(sprintf('%3.0f stimuli of the original %3.0f in the combined files will be kept in the merged output',...
            length(stims_combined_keep),length(stimulus_list)));
        if_ok=getinp('1 if ok','d',[0 1]);
    else
        stims_combined_keep=[1:length(stimulus_list)];
        if_ok=1;
    end
end %if_ok
f_base.stimulus_names_orig=strvcat(stimulus_list); %original stimulus names
%f_base.stim_labels_orig=strvcat(stim_labels_unique); %shortened names for plotting
f_base.stimulus_names=f_base.stimulus_names_orig(stims_combined_keep,:);
%f_base.stim_labels=f_base.stim_labels_orig(stims_combined_keep,:);

% Intersection is easy
resps_total = cell(2,1);
resps_total{1} = {};
resps_total{2} = {};

for setindx = 1:numSets    
    resps_concat = resps_set{setindx}(:,glomeruli_combined{1});
    resps_total{1} = [resps_total{1}; resps_concat];
end

% Union is harder
resps_tmp = resps_set;
for setindx = 1:numSets
    [numStim,numGlom] = size(resps_tmp{setindx});
    gloms_diff = setdiff(glomeruli_combined{2},resps_tmp{setindx}.Properties.VariableNames);
    resps_tmp{setindx}{:,gloms_diff} = zeros(numStim,length(gloms_diff));
    resps_concat = resps_tmp{setindx}(:,glomeruli_combined{2});
    resps_total{2} = [resps_total{2}; resps_concat];    
end


