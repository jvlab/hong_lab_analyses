function [da_new,opts_used]=hlid_da_stimselect(da,opts)
% [da_new,opts_used]=hlid_da_stimselect(da,opts) reads a raw file, typically with single-trial data,
% and keeps only a specific set of odors
%
% this enables reading of files that have validation stimuli, diagnostic stimuli, etc
% if neither opts.targets nor da specifies the list of target odors, then all stimuli are kept
%    target stimuli are obtained from da.response_amplitude_stim.is_target or da.trial_info.istarget
%    if both are present, these are checked for consistency
%
% stimuli are returned with the targets in alphabetical order
% the presence of da.trial_info is used to determine whether trial-by-trial data are present
% 
% da: a structure read from the raw file
% opts: options
%   opts.targets: cell array of specified target stimuli, e.g., {'1-5ol @ -3.0'}    {'1-6ol @ -3.0'}
%   opts.exclude_panels: panels to exclude, defaults to {'glomeruli_diagnostics'};
%   opts.if_log: 1 to log (defaults to 0)
%
% da_new: a structure with only the target stimuli kept.  The following fields are repopulated:
%   da_new.response_amplitude_stim.stim
%   da_new.response_amplitude_stim.description
%   da_new.response_amplitude_stim.mean_peak
%  and if da.trial_info is present:
%   da_new.trial_info.stim
%   da_new.response_amplitude_trials.description
%   da_new.response_amplitude_trials.baseline_win
%   da_new.response_amplitude_trials.peak_win
%   da_new.response_amplitude_trials.mean_peak
% Other fields from da are copied.
%   
%   opts_used: options used
% 
% test examples
% 
% load('C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\data\orn_validation2_wd\2023-11-19__fly01_wd.mat')
% [da_new,ou]=hlid_da_stimselect(da,setfields([],{'if_log',},{1}))
%
% da=load('C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\data\kc_tnt_singletrial\2024-07-10__fly01__TNTin_nolabel__20240710_2FLfun1.mat');
% [da_new,ou]=hlid_da_stimselect(da,setfields([],{'if_log','targets'},{1,{'1-6ol @ -3.0','ms @ -3.0'}}))
%
% See also: HLID_RASTIM2COORDS_DEMO, HLID_RASTIM2COORDS_POOL, HLID_ORN_MERGE, HLID_RASTIM_TRIAL_READ, HLID_ORN_MERGE.
%
if (nargin<=1)
    opts=struct;
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'targets',cell(0));
opts=filldefault(opts,'exclude_panels',{'glomeruli_diagnostics'});
%
% JDD, make an opt for this, so the correct field name gets returned for
% future use.
if(isfield(da.response_amplitude_stim,'mean_peak'))
    opts.numbersname='mean_peak';

elseif(isfield(da.response_amplitude_stim,'max_peak'))
    opts.numbersname='max_peak';

else
    % This could read off of a list, to which new items can be added and
    % avoid wiritng new code to add an eligible field name.
    error('Novel field name, was looking for either max_peak or mean_peak');
end

nrps = 3;

opts_used=opts;
opts_used.warnings=[];
targets_opts=opts.targets;
targets_responses=cell(0);
targets_trials=cell(0);
targets_infile=unique(da.response_amplitude_stim.stim);


% JDD 8/10/25 
% Check for rows of NaN in the raw data. Turn off the is_target flag for
% these rows.

rowsAllNaN = all(isnan(da.response_amplitude_stim.(opts.numbersname)),2);
targets_off = da.response_amplitude_stim.is_target.*rowsAllNaN';

da.response_amplitude_stim.is_target = logical(da.response_amplitude_stim.is_target - targets_off);
da.response_amplitude_stim.is_target
%prod = rowsAllNaN.*da.response_amplitude_stim.is_target

rowsAllNaN = all(isnan(da.response_amplitude_trials.(opts.numbersname)),2);
targets_off = da.trial_info.is_target.*rowsAllNaN';

da.trial_info.is_target = logical(da.trial_info.is_target-targets_off);




if isfield(da.response_amplitude_stim,'is_target')
    stim_temp=da.response_amplitude_stim.stim(da.response_amplitude_stim.is_target)
    targets_responses=unique(stim_temp);
    if length(stim_temp)>length(targets_responses)
        disp(' multiple occurrences of a target stimulus detected in stim list');
        for k=1:length(targets_responses)
            dup_ptr=strmatch(targets_responses{k},stim_temp,'exact');
            for id=2:length(dup_ptr)
                addpos=min(find(cat(2,stim_temp{dup_ptr(id)},' ')==' '));
                addstr=repmat('+',1,id-1);
                stim_temp{dup_ptr(id)}=cat(2,stim_temp{dup_ptr(id)}(1:addpos-1),addstr,stim_temp{dup_ptr(id)}(addpos:end));
            end
        end
        da.response_amplitude_stim.stim(da.response_amplitude_stim.is_target)=stim_temp;
        targets_responses=unique(stim_temp);
    end
end
if isfield(da,'trial_info') %need to detect subsequent occurrences of the same target
    if isfield(da.trial_info,'is_target')
        trial_temp=da.trial_info.stim(da.trial_info.is_target)
        targets_trials=unique(trial_temp);
        if length(trial_temp)>nrps*length(targets_trials)
            disp(' multiple occurrences of a target stimulus detected in trial list');
            for k=1:length(targets_trials)
                dup_ptr=strmatch(targets_trials{k},trial_temp,'exact');
                for id=nrps+1:length(dup_ptr)
                    addpos=min(find(cat(2,trial_temp{dup_ptr(id)},' ')==' '));
                    addstr=repmat('+',1,ceil(id/nrps)-1);
                    trial_temp{dup_ptr(id)}=cat(2,trial_temp{dup_ptr(id)}(1:addpos-1),addstr,trial_temp{dup_ptr(id)}(addpos:end));
                end
            end
            da.trial_info.stim(da.trial_info.is_target)=trial_temp;
            targets_trials=unique(trial_temp);
        end

    end
end
%
%if isfield(da.response_amplitude_stim,'is_target')
%    targets_responses=unique(da.response_amplitude_stim.stim(da.response_amplitude_stim.is_target));
%end
%if isfield(da,'trial_info')
%    if isfield(da.trial_info,'is_target')
%        targets_trials=unique(da.trial_info.stim(da.trial_info.is_target));
%    end
%end
%check that targets_responses and targets_trials match if both are present;
%if one is present, then use it
%if neither are present, use all stimuli
if_match=1;
targets_use=cell(0); %trials specified within file, or if not availble, then all stimuli in the file
if length(targets_responses)>0 & length(targets_trials)>0
    targets_use=unique([targets_responses,targets_trials]);
    if length(targets_use)>min(length(targets_responses),length(targets_trials))
        if_match=0;
    end
elseif length(targets_responses)>0
    targets_use=targets_responses;
elseif length(targets_trials)>0
    targets_use=targets_trials;
else
    targets_use=targets_infile;
end

if if_match==0
    targets_use=targets_infile;
    wstring='mismatch of target stimuli in response_amplitude_stim and trial_info, both will be ignored';
    warning(wstring);
    opts_used.warnings=strvcat(opts_used.warnings,wstring);
end
if opts.if_log | if_match==0
    disp(sprintf('unique stimuli         available  in file: %3.0f',length(targets_infile)));
    disp(sprintf('unique stimuli         specified via opts: %3.0f',length(targets_opts)));
    disp(sprintf('unique stimuli in response_amplitude_stim: %3.0f',length(targets_responses)));
    disp(sprintf('unique stimuli              in trial_info: %3.0f',length(targets_trials)));
end
if ~isempty(targets_opts)
    targets_use=targets_opts;
end
%
%for each target in target_use, determine how many occurrences in
%response_amplitude_stim.stim trial_info.stim, and in which trials
% but exclude diagnostics
resp_ptrs=zeros(1,length(targets_use));
trial_ptrs=cell(1,length(targets_use));
trial_counts=zeros(1,length(targets_use));
%
trials_exclude=[];
if isfield(da,'trial_info')
    if isfield(da.trial_info,'panel')
        for k=1:length(opts.exclude_panels)
            trials_exclude=[trials_exclude,strmatch(opts.exclude_panels{k},da.trial_info.panel,'exact')];
        end
    end
end
%
resps_exclude=[];
if isfield(da.response_amplitude_stim,'panel')
    for k=1:length(opts.exclude_panels)
        resps_exclude=[resps_exclude,strmatch(opts.exclude_panels{k},da.response_amplitude_stim.panel,'exact')];
    end
end
for istim=1:length(targets_use)
    resp_ptr=strmatch(targets_use{istim},da.response_amplitude_stim.stim,'exact');
    resp_ptr=setdiff(resp_ptr,resps_exclude);
    if length(resp_ptr)~=1
        wstring=sprintf('requested stimulus %s not found or multiply found in response_amplitude_stim.stim',targets_use{istim});
        warning(wstring);
        opts_used.warnings=strvcat(opts_used.warnings,wstring);
    else
        resp_ptrs(istim)=resp_ptr;
    end
    if isfield(da,'trial_info')
        trial_ptr=strmatch(targets_use{istim},da.trial_info.stim,'exact');
        trial_ptr=setdiff(trial_ptr,trials_exclude);
         if length(trial_ptr)==0
            wstring=sprintf('requested stimulus %s not found in trial_info.stim',targets_use{istim});
         else
            trial_ptrs{istim}=trial_ptr;
            trial_counts(istim)=length(trial_ptr);
        end
    end
end
opts_used.resp_ptrs=resp_ptrs;
opts_used.trial_ptrs=trial_ptrs;
opts_used.trial_counts=trial_counts;
%consistency checks done, pull requested stimuli and trials
da_new=da;
%response amplitudes
ras=struct();
ras.description=da.response_amplitude_stim.description;
ras.(opts.numbersname)=NaN(length(targets_use),size(da.response_amplitude_stim.(opts.numbersname),2));
ras.stim=cell(1,length(targets_use));
for istim=1:length(targets_use)
    if resp_ptrs(istim)>0
        ras.(opts.numbersname)(istim,:)=da.response_amplitude_stim.(opts.numbersname)(resp_ptrs(istim),:);
    end
    ras.stim{istim}=targets_use{istim};
end
da_new.response_amplitude_stim=ras;
%trial-by-trial data
if isfield(da,'trial_info')
    rat=struct;
    rat.description=da.response_amplitude_trials.description;
    rat.baseline_win=da.response_amplitude_trials.baseline_win;
    rat.peak_win=da.response_amplitude_trials.peak_win;
    rat.(opts.numbersname)=NaN(sum(trial_counts),size(da.response_amplitude_trials.(opts.numbersname),2));
    ti=struct;
    ti.stim=cell(1,sum(trial_counts));
    itrial=0;
    for istim=1:length(targets_use)
        for k=1:trial_counts(istim)
            itrial=itrial+1;
            rat.(opts.numbersname)(itrial,:)=da.response_amplitude_trials.(opts.numbersname)(trial_ptrs{istim}(k),:);
            ti.stim{itrial}=targets_use{istim};
        end
    end
    da_new.response_amplitude_trials=rat;
    da_new.trial_info=ti;
end
return
end
