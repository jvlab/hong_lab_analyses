%hlid_pred_magnif_demo: predict distances in transformation from orn to kc
%
% reads a set of orn files, checks for consistency, creates a merged file from trial-averaged z-scores
% via filling in with multiplicative and additive offset (via afalwt, as in hlid_orn_merge, with if_restore_size=1)
%
%   See also:  HLID_SETUP, HLID_ORN_MERGE, HLID_FILL_MERGE_SVD, HLID_COORDS_SVD.
%
hlid_setup;
if ~exist('opts_dasel')
    opts_dasel=struct;
end
if_restore_size=getinp('1 to restore responses to match overall size of original data (0: legacy)','d',[0 1],1);
if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
%
[filenames_short,pathname]=uigetfile('*fly*.mat','Select raw ORN data files','Multiselect','on');
if ~iscell(filenames_short)
    filenames_short=cellstr(filenames_short);
end
nfiles=length(filenames_short);
s=cell(nfiles,1);
files_use=[];
nancols=cell(nfiles,1);
nanrows=cell(nfiles,1);
dsid=[];
%
%verify consistency of names and numbers of glomeruli and stimuli
%
for ifile=1:nfiles
    s{ifile}=load(cat(2,pathname,filenames_short{ifile}));
    [s{ifile},optsused_dasel]=hlid_da_stimselect(s{ifile},opts_dasel);
    dsid_this=s{ifile}.meta.title;
    %'-'can be used within fields of file name
    dsid_this=strrep(dsid_this,'/','-');
    dsid_this=strrep(dsid_this,'\','-');
    dsid_this=strrep(dsid_this,'_','-');
    while contains(dsid_this,'--')
        dsid_this=strrep(dsid_this,'--','-');
    end
    %
    dsid=strvcat(dsid,dsid_this);
    if (ifile==1)
        glomeruli=s{ifile}.rois.glomeruli;
        nglomeruli=length(glomeruli);
        stimulus_names=s{ifile}.response_amplitude_stim.stim';
        nstims=length(stimulus_names);
    end
    glomeruli_check=s{ifile}.rois.glomeruli;
    stimulus_names_check=s{ifile}.response_amplitude_stim.stim';
    resps_raw=s{ifile}.response_amplitude_stim.mean_peak;
    nancols{ifile}=find(all(isnan(resps_raw),1));
    nanrows{ifile}=find(all(isnan(resps_raw),2));
    %
    disp(sprintf('file %2.0f (%20s) read, %3.0f glomeruli (%3.0f all NaN), %3.0f stimuli (%3.0f all NaN)',ifile,filenames_short{ifile},...
        length(glomeruli_check),length(nancols{ifile}),...
        length(stimulus_names_check),length(nanrows{ifile})));
    ifok=1;
    if length(glomeruli_check)~=nglomeruli
        ifok=0;
        disp('number of glomeruli does not match')
    elseif any(strcmp(glomeruli,glomeruli_check)==0)
        ifok=0;
        disp('names of glomeruli do not match')
    end
    if length(stimulus_names_check)~=nstims
        ifok=0;
        disp('number of stimuli does not match')
    elseif any(strcmp(stimulus_names,stimulus_names_check)==0)
        ifok=0;
        disp('names of stimuli do not match')
    end
    if (ifok==1)
        files_use=[files_use,ifile];
    end
end
%shorten stimulus names
stim_labels=stimulus_names;
for istim=1:nstims
    if contains(stim_labels{istim},'@')
        stim_labels{istim}=stim_labels{istim}(1:min(find(stim_labels{istim}=='@')-1));
    end
    stim_labels{istim}=deblank(stim_labels{istim});
end
%
% select datasets
%
files_use=getinp('list of files to use','d',[1 nfiles],files_use);
nfiles_use=length(files_use);
min_present=getinp('minimum number of preps that a glomerulus must be present in','d',[0 nfiles_use],1);
%
opts_fill_merge=struct;
opts_fill_merge.files_use=files_use;
opts_fill_merge.glomeruli=glomeruli;
opts_fill_merge.nancols=nancols; %this needs to be determined on the basis of all stimuli
opts_fill_merge.min_present=min_present;
opts_fill_merge.if_log=1;
opts_fill_merge.if_restore_size=if_restore_size;
opts_fill_merge.stim_labels=stim_labels;
opts_fill_merge.response_name='response_amplitude_stim';
opts_fill_merge.response_name2='mean_peak';
opts_fill_merge.stimuli_use=[1:nstims];
opts_fill_merge.if_plot=1;
opts_fill_merge.filenames_short=filenames_short;
opts_fill_merge.if_restore_size=if_restore_size;
opts_fill_merge.if_submean=if_submean;
%
[resps,coords_all]=hlid_fill_merge_svd(s,opts_fill_merge);
if any(isnan(resps(:)))
    disp('Cannot proceed. Not all NaNs have been filled in.')
end
%
%import the data into rs format
%
nds=size(coords_all,2);
coords=cell(1,nds);
for id=1:nds
    coords{id}=coords_all(:,1:id);
end
aux_import=struct;
opts_import=struct;
opts_import.typenames=stim_labels;
opts_import.type_coords_def='none';
opts_import.paradigm_name='orn terminals';
opts_import.label_long=cat(2,'orn terminals trial-averaged, merged, svd',sprintf(' meansub=%1.0f',if_submean));
aux_import.opts_import=opts_import;
[data_orn,aux_import_out]=rs_import_coordsets(coords,aux_import);


%this is to test that drop-2 works%
opts_fill_merge_drop2=opts_fill_merge;
opts_fill_merge_drop2.stimuli_use=[2:16];
[resps_drop2,coords_all_drop2]=hlid_fill_merge_svd(s,opts_fill_merge_drop2);
%

%also need to import a KC dataset
