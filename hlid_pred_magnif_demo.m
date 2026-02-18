%hlid_pred_magnif_demo: predict distances in transformation from orn to kc
%
% ORN: reads a set of raw data, trial-averaged files, checks for consistency, creates a merged file from trial-averaged z-scores
% via filling in with multiplicative and additive offset (via afalwt, as in hlid_orn_merge, with if_restore_size=1).
% The missing-data procedure is carried out before any stimuli are dropped.
%
% KC: reads a set of raw data, trial-averaged files, for each, creates a coordinate set via svd, and combines
% via Procrustes consensus with offset and without scaling
%
% It is necessary to start from raw data files since the prediction is based on leaving out two stimuli,
% so the spaces need to be constructed from scratch with those stimuli deleted
%
% mean subtraction (optional, determined by if_submean, defaults to 0):
%   If mean subtraction is enabled, for the ORN data, it is carried out after the preps are merged by the missing-data SVD routine
%   For KC data, it is carried out within each prep, prior to svd, and prior to Procrustes
%
% For fitting the transformation from ORN space to KC space, the ORN coordinaates are restricted to the retained stmuli:
%   coordinates are coords_all_orn, computed from resps_orn (responses after dropping stimuli, merging by missing-data method) in hlid_merge_svd.
%   These are re-computed by coords_orn_drop=resps_orn*v_drop, where v_drop is the projection from ORN responses to the svd space
%   and we verify that coords_orn_drop=coords_all_orn (dev_orn_coords, dev_max_orn_coords)
% For predicting magnification factors, coords_orn_all_drop is the stimuli, mapped into the the dropped-stimulus space,
%  via the same projection v_drop: coords_orn_all_drop=resps_orn_full*v_drop, and the affine transformation is applied to it.
%  Note that the coordinates of the *un*dropped stimuli in resps_orn_full are not identical to the responses of the same stimuli in resps_orn,
%   because of the way that amplitudes are matched in hlid_fill_merge_svd, but this difference is small (dev_orn_coords_drop, dev_max_orn_coords_drop)
%
%   Alternative (not done) would be to use just the merge of all stimuli, but just the coordinates of the un-dropped stimuli.
%   This could be done by selecting the coordinates from coords_all_orn_full, which is coords_all_orn as computed without any dropped stimuli.
%
%   See also:  HLID_SETUP, HLID_ORN_MERGE, HLID_FILL_MERGE_SVD, HLID_COORDS_SVD, HLID_DA_STIMSELECT, RS_GEOFIT.
%
hlid_setup;
if ~exist('opts_dasel')
    opts_dasel=struct;
end
if_restore_size=getinp('1 to restore responses to match overall size of original data (0: legacy)','d',[0 1],1);
if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
%
[filenames_short_orn,pathname_orn]=uigetfile('*fly*.mat','Select raw ORN data files','Multiselect','on');
if ~iscell(filenames_short_orn)
    filenames_short_orn=cellstr(filenames_short_orn);
end
nfiles_orn=length(filenames_short_orn);
s_orn=cell(nfiles_orn,1);
files_use_orn=[];
nancols_orn=cell(nfiles_orn,1);
nanrows_orn=cell(nfiles_orn,1);
dsid=[];
%
%verify consistency of names and numbers of glomeruli and stimuli
%
for ifile=1:nfiles_orn
    s_orn{ifile}=load(cat(2,pathname_orn,filenames_short_orn{ifile}));
    [s_orn{ifile},optsused_dasel]=hlid_da_stimselect(s_orn{ifile},opts_dasel);
    dsid_this=s_orn{ifile}.meta.title;
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
        glomeruli=s_orn{ifile}.rois.glomeruli;
        nglomeruli=length(glomeruli);
        stimulus_names=s_orn{ifile}.response_amplitude_stim.stim';
        nstims=length(stimulus_names);
    end
    glomeruli_check=s_orn{ifile}.rois.glomeruli;
    stimulus_names_check=s_orn{ifile}.response_amplitude_stim.stim';
    resps_raw=s_orn{ifile}.response_amplitude_stim.mean_peak;
    nancols_orn{ifile}=find(all(isnan(resps_raw),1));
    nanrows_orn{ifile}=find(all(isnan(resps_raw),2));
    %
    disp(sprintf('file %2.0f (%20s) read, %3.0f glomeruli (%3.0f all NaN), %3.0f stimuli (%3.0f all NaN)',ifile,filenames_short_orn{ifile},...
        length(glomeruli_check),length(nancols_orn{ifile}),...
        length(stimulus_names_check),length(nanrows_orn{ifile})));
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
        files_use_orn=[files_use_orn,ifile];
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
files_use_orn=getinp('list of files to use for ORN data','d',[1 nfiles_orn],files_use_orn);
nfiles_use_orn=length(files_use_orn);
min_present=getinp('minimum number of preps that a glomerulus must be present in','d',[0 nfiles_use_orn],1);
%
%Kenyon cell  data
%
[filenames_short_kc,pathname_kc]=uigetfile('*fly*megamat*.mat','Select raw KC data files','Multiselect','on');
if ~iscell(filenames_short_kc)
    filenames_short_kc=cellstr(filenames_short_kc);
end
nfiles_kc=length(filenames_short_kc);
s_kc=cell(nfiles_kc,1);
files_use_kc=[];
nancols_kc=cell(nfiles_kc,1);
nanrows_kc=cell(nfiles_kc,1);
dsid_kc=[];
%
%verify consistency of names with ORN files
%
for ifile=1:nfiles_kc
    s_kc{ifile}=load(cat(2,pathname_kc,filenames_short_kc{ifile}));
    [s_kc{ifile},optsused_dasel_kc]=hlid_da_stimselect(s_kc{ifile},opts_dasel);
    dsid_this=s_kc{ifile}.meta.title;
    %'-'can be used within fields of file name
    dsid_this=strrep(dsid_this,'/','-');
    dsid_this=strrep(dsid_this,'\','-');
    dsid_this=strrep(dsid_this,'_','-');
    while contains(dsid_this,'--')
        dsid_this=strrep(dsid_this,'--','-');
    end
    %
    dsid_kc=strvcat(dsid_kc,dsid_this);
    stimulus_names_check=s_kc{ifile}.response_amplitude_stim.stim';
    stimulus_names_check=strrep(stimulus_names_check,'.0',''); %KC data has concentration labels 3.0, not 3
    resps_raw=s_kc{ifile}.response_amplitude_stim.mean_peak;
    nancols_kc{ifile}=find(all(isnan(resps_raw),1));
    nanrows_kc{ifile}=find(all(isnan(resps_raw),2));
    disp(sprintf('file %2.0f (%20s) read, %3.0f rois (%3.0f all NaN), %3.0f stimuli (%3.0f all NaN)',ifile,filenames_short_orn{ifile},...
        size(resps_raw,2),length(nancols_kc{ifile}),...
        length(stimulus_names_check),length(nanrows_kc{ifile})));
    ifok=1;
    if length(stimulus_names_check)~=nstims
        ifok=0;
        disp('number of stimuli does not match')
    elseif any(strcmp(stimulus_names,stimulus_names_check)==0)
        ifok=0;
        disp('names of stimuli do not match')
    end
    if (ifok==1)
        files_use_kc=[files_use_kc,ifile];
    end
end
files_use_kc=getinp('list of files to use for KC data','d',[1 nfiles_kc],files_use_kc);
nfiles_use_kc=length(files_use_kc);
%
dim_max=getinp('maximum analysis dimension','d',[1 nstims],7);
drop_stim_max=getinp('maximum stimulus number to drop','d',[2 nstims],nstims);  %for debugging
%
opts_fill_merge=struct;
opts_fill_merge.files_use=files_use_orn;
opts_fill_merge.glomeruli=glomeruli;
opts_fill_merge.nancols=nancols_orn; %this needs to be determined on the basis of all stimuli
opts_fill_merge.min_present=min_present;
opts_fill_merge.if_log=1;
opts_fill_merge.if_restore_size=if_restore_size;
opts_fill_merge.stim_labels=stim_labels;
opts_fill_merge.response_name='response_amplitude_stim';
opts_fill_merge.response_name2='mean_peak';
opts_fill_merge.if_plot=1;
opts_fill_merge.filenames_short=filenames_short_orn;
opts_fill_merge.if_restore_size=if_restore_size;
opts_fill_merge.if_submean=if_submean;
%
opts_import_orn_all=struct;
opts_import_orn_all.paradigm_name='orn terminals';
opts_import_orn_all.label_long=cat(2,'orn terminals trial-averaged, merged, svd',sprintf(' meansub=%1.0f',if_submean));
opts_import_orn_all.label='orn';
opts_import_orn_all.type_coords_def='none';
%
opts_import_kc_all=struct;
opts_import_kc_all.type_coords_def='none';
opts_import_kc_all.paradigm_name='kc soma';
%
% loop over stims_drop=[] and all pairs
% if [], then save the merged orn and kc files separately
% if not, then log that the pair is done
%
% in each case,
% project each glomerular pattern onto the coordinates, taking into account if_submean
% do this with the glomerular pattern from the full set (which includes the stimuli left in), and also the left out stimuli
% the left-in stimuli simply serve as a check
%
% then fit the affine model, and also, need to save how
% the glomeruli are projected into the coord space -- including taking into account if_submean
%
drop_list=nchoosek([1:drop_stim_max],2);
ndrop_list=size(drop_list,1);
dev_max_orn_coords=0;
dev_max_orn_coords_drop=0;
d_fits=zeros(ndrop_list+1,dim_max); %goodness of fits
for idrop=0:ndrop_list
    if idrop==0
        stims_drop=[];
        drop_text='full set';
    else
        stims_drop=drop_list(idrop,:);
        drop_text=sprintf('set with dropped stimuli [%2.0f %2.0f]',stims_drop);
    end
    stims_use=setdiff([1:nstims],stims_drop);
    %
    %merge and import the ORN data into rs format
    %
    opts_fill_merge.stimuli_use=stims_use;
    opts_fill_merge.if_log=double(isempty(stims_drop));
    opts_fill_merge.if_plot=double(isempty(stims_drop));
    [resps_orn,coords_all_orn,f_orn,resps_mean_orn]=hlid_fill_merge_svd(s_orn,opts_fill_merge);
    if_notok=any(isnan(resps_orn(:)));
    if if_notok
        disp(sprintf('skipping %s, not all NaNs have been filled in.',drop_text));
        if isempty(stims_dropped)
            disp('Cannot proceed');
        end
    else
        %
        opts_import_orn_all.typenames=stim_labels(stims_use);
        aux_import=struct;
        aux_import.opts_import=opts_import_orn_all;
        [data_orn_all,aux_import_orn_all]=rs_import_coordsets(coords_all_orn,aux_import);
        %
        %for KC apply PCA to get coordinates, and import each set separately
        %
        data_kc_all=cell(1,nfiles_use_kc);
        for ifile_ptr=1:nfiles_use_kc
            ifile=files_use_kc(ifile_ptr);
            resps_raw=s_kc{ifile}.response_amplitude_stim.mean_peak;
            stims_keep=setdiff([1:nstims],union(nanrows_kc{ifile},stims_drop));
            rois_keep=setdiff([1:size(resps_raw,2)],nancols_kc{ifile});
            resps_kc=s_kc{ifile}.response_amplitude_stim.mean_peak(stims_keep,rois_keep);           
            if if_submean %subtract the mean of the stimuli being analyzed               
                resps_kc=resps_kc-repmat(mean(resps_kc,1),size(resps_kc,1),1);
            end
            maxdim_allowed=min(size(resps_kc))-if_submean;
            maxdim_use=maxdim_allowed;
            [f,s_diag_all,u_full,v_full,s_full,coords_all_kc]=hlid_coords_svd(struct(),resps_kc,maxdim_allowed,maxdim_use,if_submean,[],[],setfield(struct(),'if_log',0));
            %
            opts_import_kc_all.typenames=stim_labels(stims_keep);
            opts_import_kc_all.label_long=cat(2,'kc soma trial-averaged, svd ',filenames_short_kc{ifile},sprintf(' meansub=%1.0f',if_submean));
            opts_import_kc_all.label=cat(2,'kc ', strrep(filenames_short_kc{ifile},'.mat',''));
            aux_import=struct;
            aux_import.opts_import=opts_import_kc_all;
            [data_kc_onefile,aux_import_kc_all]=rs_import_coordsets(coords_all_kc,aux_import);
            if ifile_ptr==1
                data_kc_all=data_kc_onefile;
            else
                data_kc_all=rs_concat_coordsets(data_kc_all,data_kc_onefile);
            end
            aux_knit=struct;
            aux_knit.opts_knit.dim_max_in=dim_max;
            aux_knit.opts_knit.if_log=double(isempty(stims_drop));
            aux_knit.opts_check.if_warn=double(isempty(stims_drop));
            [data_kc_knit,aux_knit_out]=rs_knit_coordsets(data_kc_all,aux_knit);
        end
        if idrop==0
            drawnow;
            %save values from analysis with full stimulus set
            resps_orn_full=resps_orn;
            resps_mean_orn_full=resps_mean_orn;
            coords_all_orn_full=coords_all_orn;
            svd_orn_full=f_orn.coord_opts;
            data_orn_full=data_orn_all;
            data_kc_full=data_kc_knit;
        end
        svd_orn=f_orn.coord_opts;
        %
        %determine affine model from dropped ORN coords to dropped KC coords
        %
        aux_geof=struct;
        aux_geof.opts_geof.model_list={'affine_offset'};
        aux_geof.opts_geof.if_fit_summary=double(isempty(stims_drop));
        aux_geof.opts_geof.if_warn=double(isempty(stims_drop));
        aux_geof.opts_geof.if_log=double(isempty(stims_drop));
        [gfs,xs,aux_geof_out]=rs_geofit(data_orn_all,data_kc_knit,aux_geof);
        for idim=1:dim_max
            d_fits(idrop+1,idim)=gfs{1}.gf{idim,idim}.d;
        end
        %
        %compute coordinates in ORN space of all stimuli,including dropped ones
        %
        v_full=svd_orn_full.aux.v(:,1:dim_max);
        v_drop=f_orn.coord_opts.aux.v(:,1:dim_max);
        coords_orn_all_full=resps_orn_full*v_full; %coordinates obtained from merging all stimuli, mapped into ORN dataset with all stimuli
        coords_orn_all_drop=resps_orn_full*v_drop; %coordinates obtained from merging all stimuli, mapped into ORN dataset with dropped stimuli
        coords_orn_drop=resps_orn*v_drop; %coordinates obtained from merging stimuli after drop, mapped into ORN dataset with dropped stimuli, should match coords_all_orn;
        dev_orn_coords=max(max(abs(coords_orn_drop-coords_all_orn(:,1:dim_max))));
        dev_max_orn_coords=max(dev_max_orn_coords,dev_orn_coords);
        dev_orn_coords_drop=max(max(abs(coords_orn_all_drop(stims_use,:)-coords_orn_drop)));
        dev_max_orn_coords_drop=max(dev_max_orn_coords_drop,dev_orn_coords_drop);
        disp(sprintf('analyzed %s, orn cooord consistency check: %10.8f, drop effect: %10.8f',drop_text,dev_orn_coords,dev_orn_coords_drop));
     end %if_notok
end %drop_list
disp(sprintf('max orn cooord consistency check: %10.8f, drop effect: %10.8f',dev_max_orn_coords,dev_max_orn_coords_drop));
disp('goodness of fits for affine models for full dataset, as function of dimension');
disp(d_fits(1,:));
disp('mean goodness of fits for affine models for datasets with dropped stimuli');
disp(mean(d_fits(2:end,:),1));

