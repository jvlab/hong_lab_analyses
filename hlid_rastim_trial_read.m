%hlid_rastim_trial_read script to red read calcium imaging data from Hong Lab at single-trial level
%
% 23Feb255: add if_log, can set to 0 to suppress most output
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM_TRIAL_PCA, HLID_RASTIM_TRIAL_VIS.
% 
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/data/kc_soma_nls/2022-10-10__fly01__megamat0.mat'; end
if ~exist('if_log') if_log=1; end
nfiles_signed=0;
while nfiles_signed==0
    nfiles_signed=getinp('number of files to analyze (<0 for dialog box)','d',[-Inf Inf]);
end
nfiles=abs(nfiles_signed);
%
metadata=cell(nfiles,1);
das=cell(nfiles,1);
dsids=cell(nfiles,1);
resps_mean=cell(nfiles,1); %mean responses across stimuli
trial_ptrs=cell(nfiles,1); %array of (nstims,nrepts)
resps_trial=cell(nfiles,1); %trial-by-trial responses, stimuli unscrambled, (nstims, nrepts, nrois_avail)
trial_sequence=cell(nfiles,1); %stimulus sequence, stimuli as strings
stims_avail=cell(nfiles,1); %list of available stimuli in each file, beginning at 1
rois_avail=cell(nfiles,1); %list of roi numbers kept for analysis, beginning at 1
rois=cell(nfiles,1); % original rois
%
nrois_avail=zeros(1,nfiles);
%
if_ok=0;
ui_prompt='Select Hong lab files to be pooled';
ui_filter_base='_fly';
if ~exist('raw_path') raw_path='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/data/'; end
ui_filter=cat(2,raw_path,'*',ui_filter_base,'*.mat');
while (if_ok==0)
    if (nfiles_signed<0)
        [filenames_short,pathname,filter_index]=uigetfile(ui_filter,ui_prompt,'Multiselect','on');
    end
    %non-UI method
    for ifile=1:nfiles
        if (nfiles_signed>0)
            HongLab_fn=getinp(sprintf('Hong Lab file name for file %2.0f to be analyzed',ifile),'s',[],HongLab_fn);
        else
            if ~iscell(filenames_short)
                filenames_short=cellstr(filenames_short);
            end
            disp(sprintf(' file %2.0f: %s',ifile,filenames_short{ifile}));
            HongLab_fn=cat(2,pathname,filesep,filenames_short{ifile});
        end
        das{ifile}=load(HongLab_fn);
        metadata{ifile}=das{ifile}.meta;
        if isfield(das{ifile},'rois')
            rois{ifile}=das{ifile}.rois;
        end
        stimulus_names_this=strvcat(das{ifile}.response_amplitude_stim.stim');
        dsid_this=das{ifile}.meta.title;
        %'-'can be used within fields of file name
        dsid_this=strrep(dsid_this,'/','-');
        dsid_this=strrep(dsid_this,'\','-');
        dsid_this=strrep(dsid_this,'_','-');
        while contains(dsid_this,'--')
            dsid_this=strrep(dsid_this,'--','-');
        end
        %
        dsids{ifile}=dsid_this;
        nstims_this=length(stimulus_names_this);
        if ifile==1
            nstims=nstims_this;
            stimulus_names=stimulus_names_this;
        end
        if_ok=1;
        if nstims~=nstims_this %check that the stimulus list is identical
            if_ok=0;
        else
            if any(stimulus_names~=stimulus_names_this)
                if_ok=0;
            end
        end
        %
        if (if_ok==0)
            disp('incompatible number of stimuli or stimulus names');
        else %file has compatible stimuli; now check dimensions
            resps_mean{ifile}=das{ifile}.response_amplitude_stim.mean_peak;
            if size(resps_mean{ifile},1)~=nstims
                warning(sprintf('number of responses (%2.0f) is not equal to number of stimuli (%2.0f)',size(resps,1),nstims));
            end
            %
            trial_sequence{ifile}=das{ifile}.trial_info.stim;
              if length(trial_sequence{ifile})~=nrepts*nstims
                warning(sprintf('number of trials (%2.0f) is not number of stims (%2.0f) * number of expected repeats (%2.0f)',...
                    length(trial_sequence{ifile}),nstims,nrepts));
            end
            %are all stimuli found?
            stim_ptrs=zeros(1,length(trial_sequence{ifile}));
            for itrial=1:length(trial_sequence{ifile})
                istim=strmatch(trial_sequence{ifile}{itrial},stimulus_names,'exact');
                if isempty(istim)
                    istim=strmatch(strrep(trial_sequence{ifile}{itrial},'.0',''),stimulus_names,'exact');
                end
                if isempty(istim)
                    warning(sprintf('trial %2.0f: stimulus label %s not recognized',itrial,trial_sequence{ifile}{itrial}));
                else
                    stim_ptrs(itrial)=istim;
                end
            end
            %are all stimuli present in nrepts trials?
            trial_ptrs{ifile}=zeros(nstims,nrepts);
            for istim=1:nstims
                trials_found=find(stim_ptrs==istim);
                if length(trials_found)~=nrepts
                    warning(sprintf('stimulus %2.0f found in %2.0f trials, expecting %2.0f',istim,length(trials_found),nrepts));
                else
                    trial_ptrs{ifile}(istim,:)=trials_found;
                end
                if if_log
                    disp(cat(2,sprintf(' stimulus %2.0f (%s) appears in trials',istim,stimulus_names(istim,:)),sprintf(' %4.0f',trials_found)));
                end
            end
            %
            nancols_all=find(all(isnan(resps_mean{ifile}),1));
            nancols_any=find(any(isnan(resps_mean{ifile}),1));
            nanrows_all=find(all(isnan(resps_mean{ifile}),2));
            stims_avail{ifile}=setdiff([1:nstims],nanrows_all);
            nancols_sel=find(any(isnan(resps_mean{ifile}(stims_avail{ifile},:)),1));
            disp(sprintf('ROIs found: %3.0f, ROIs that are all NaNs: %3.0f; ROIs that have any NaNs: %3.0f; ROIs with NaNs when stimuli are present: %3.0f',...
                size(resps_mean{ifile},2),length(nancols_all),length(nancols_any),length(nancols_sel)));
            disp(sprintf('unique stimuli found: %4.0f',nstims));
            rois_avail{ifile}=setdiff([1:size(resps_mean{ifile},2)],nancols_sel);
            nrois_avail(ifile)=length(rois_avail{ifile});
            resps_mean{ifile}=resps_mean{ifile}(:,rois_avail{ifile}); %remove ROI's with missing data but retain stimuli with all-NaN responses
          end
    end %ifile
    if (if_ok)
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end %if_ok
%sort trials according to stimuli and summarize
for ifile=1:nfiles
    for istim=1:nstims
        resps_trial{ifile}(istim,:,:)=reshape(das{ifile}.response_amplitude_trials.mean_peak(trial_ptrs{ifile}(istim,:),rois_avail{ifile}),[1 nrepts nrois_avail(ifile)]);
    end
    check_mean=max(max(abs(resps_mean{ifile}-squeeze(mean(resps_trial{ifile},2)))));
    disp(sprintf(' file %2.0f (%30s):  mean responses: %4.0f x %4.0f, stims available: %4.0f, compare computed trial mean and supplied mean: %9.7f',...
        ifile,dsids{ifile},size(resps_mean{ifile}),length(stims_avail{ifile}),check_mean));
end
clear das

