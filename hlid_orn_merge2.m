%hlid_orn_merge2: read multiple sets of raw orn data files with non-overlapping stimuli
% This handles missing data, merges the orn reponses, and writes a coordinate file.
%
% Note that its inputs are raw data files, not coordinate files -- so it
% does not use the psg library for file input, and cannot use psg_align_coordsets to align.
%
% Missing data within a set are filled in as in hlid_orn_merge,
% but then the sets are combined, using corresponding ORNs.
% 
% If the sets are non-overlapping, they are appended.  If there are responses in common, they are averaged.
%
% rois.glomeruli fields must match across sets
%
% Built on hlid_orn_merge
% 
% * many quantities of hlid_orn_merge are here arrays or cell arrays of size (1,nsets)
% * each set of files must have non-overlapping stimuli
% * merges across sets of files (sets of odor panels) by either an intersection (glomeruli present in all) or
%     union (glomeruli present in any,, with missing glomeruli set to zero
% * defaults to restoring responses to match overall size of original data (hlid_orn_merge does not)
% * merged coordinate sets always go up to largest possible dimension (hlid_orn_merge allows choice)
%
% 03Jun25: allow for the sets to have overlapping stimuli. Logic changed: responses are added in rather than concatenated.
% 01Jul25: add roi_names and glomerulus labels to outputs
%
%   See also:  HLID_SETUP, HLID_RASTM2COORDS_DEMO, HLID_RASTIM2COORDS_POOL, AFALWT, HLID_SVD_COORDS, HLID_PLOT_COORDS,
%   HLID_DA_STIMSELECT, HLID_ORN_MERGE.
%
hlid_setup;
if ~exist('opts_dasel')
    opts_dasel=struct;
end
if ~exist('hist_bins')  hist_bins=50; end
if ~exist('hist_quantiles') hist_quantiles=[.05 .25 .5 .75 .95]; end
nsets=getinp('number of sets of stimuli','d',[1 100]);
if_restore_size=getinp('1 to restore responses to match overall size of original data (0: legacy)','d',[0 1],1);
if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
%
filenames_short=cell(1,nsets);
pathname_short=cell(1,nsets);
nfiles=zeros(1,nsets);
nflies_use=zeros(1,nsets);
s=cell(1,nsets);
stimulus_names_set=cell(1,nsets);
stim_labels_set=cell(1,nsets);
nstims=zeros(1,nsets);
glomeruli_use=cell(1,nsets); %pointers to glomeruli used in each set
nglomeruli_use=zeros(1,nsets); %number of glomeruli used in each set
resps_raw=cell(1,nsets); %responses within each dataset prior to fill-in merge
files_use=cell(1,nsets); %files used within each set
dsid=cell(1,nsets);
%
for iset=1:nsets
    disp(sprintf('Enter set %1.0f',iset))
    [filenames_short{iset},pathname{iset}]=uigetfile('*fly*.mat',sprintf('Select raw ORN data files for set %1.0f',iset),'Multiselect','on');
    if ~iscell(filenames_short{iset}) % Why is this necessary? When would this happen? Does it happen every time? Does it matter?
        filenames_short{iset}=cellstr(filenames_short{iset});
    end
    nfiles(iset)=length(filenames_short{iset}); % The number of files that comprise this set
    files_use{iset}=[]; % Don't know what this is for yet.
    nancols=cell(nfiles(iset),1); % I can guess, but I won't
    nanrows=cell(nfiles(iset),1);
    dsid{iset}=[]; % I am guessing (d)ata (s)et (id)entification
    %
    %verify consistency of names and numbers of glomeruli and stimuli
    %
    for ifile=1:nfiles(iset)
        s{iset}{ifile}=load(cat(2,pathname{iset},filenames_short{iset}{ifile})); % Loads the data from the file. S holds the familiar data structure.
        [s{iset}{ifile},optsused_dasel]=hlid_da_stimselect(s{iset}{ifile},opts_dasel);
        dsid_this=s{iset}{ifile}.meta.title;
        %'-'can be used within fields of file name
        dsid_this=strrep(dsid_this,'/','-');
        dsid_this=strrep(dsid_this,'\','-');
        dsid_this=strrep(dsid_this,'_','-');
        while contains(dsid_this,'--')
            dsid_this=strrep(dsid_this,'--','-');
        end
        %
        dsid{iset}=strvcat(dsid{iset},dsid_this);
        if (ifile==1)     % Why is this here? One of two things needs to be done.
            if (iset==1)  % Either lose this requirement - why do we care if the sets are consistent?
                glomeruli=s{iset}{ifile}.rois.glomeruli;
                nglomeruli=length(glomeruli);
            end
            stimulus_names_set{iset}=s{iset}{ifile}.response_amplitude_stim.stim';
            nstims(iset)=length(stimulus_names_set{iset});
        end
        glomeruli_check=s{iset}{ifile}.rois.glomeruli;
        stimulus_names_check=s{iset}{ifile}.response_amplitude_stim.stim';
        resps_read=s{iset}{ifile}.response_amplitude_stim.(optsused_dasel.numbersname);
        nancols{ifile}=find(all(isnan(resps_read),1));
        nanrows{ifile}=find(all(isnan(resps_read),2));% The NaN situation is addressed here, though there doesn't seem to be anything done about it below.
        %
        disp(sprintf('file %2.0f (%20s) read, %3.0f glomeruli (%3.0f all NaN), %3.0f stimuli (%3.0f all NaN)',ifile,filenames_short{iset}{ifile},...
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
        if length(stimulus_names_check)~=nstims(iset)
            ifok=0;
            disp('number of stimuli does not match')
        elseif any(strcmp(stimulus_names_set{iset},stimulus_names_check)==0)
            ifok=0;
            disp('names of stimuli do not match')
        end
        if (ifok==1)
            files_use{iset}=[files_use{iset},ifile];
        end
    end
    %shorten stimulus names
    stim_labels_set{iset}=stimulus_names_set{iset};
    for istim=1:nstims(iset)
        if contains(stim_labels_set{iset}{istim},'@')
            stim_labels_set{iset}{istim}=stim_labels_set{iset}{istim}(1:min(find(stim_labels_set{iset}{istim}=='@')-1));
        end
        stim_labels_set{iset}{istim}=deblank(stim_labels_set{iset}{istim});
    end
    %
    % select datasets
    %
    files_use{iset}=getinp('list of files to use','d',[1 nfiles(iset)],files_use{iset});
    nfiles_use(iset)=length(files_use{iset});
    resps_raw{iset}=zeros(nstims(iset),nglomeruli,nfiles_use(iset));
    glomeruli_missing=zeros(nglomeruli,nfiles_use(iset));
    for ifile_ptr=1:nfiles_use(iset)
        ifile=files_use{iset}(ifile_ptr);
        glomeruli_missing(nancols{ifile},ifile_ptr)=1;
        resps_raw{iset}(:,:,ifile_ptr)=s{iset}{ifile}.response_amplitude_stim.(optsused_dasel.numbersname);
    end
    for ipresent=0:nfiles_use(iset)
        disp(sprintf('number of glomeruli present in at least %2.0f datasets: %3.0f (missing in %2.0f or less)',...
            ipresent,sum(sum(glomeruli_missing,2)<=(nfiles_use(iset)-ipresent)),nfiles_use(iset)-ipresent));
    end
    min_present=getinp('minimum number of preps that a glomerulus must be present in','d',[0 nfiles_use(iset)],1);
    glomeruli_use{iset}=find(sum(glomeruli_missing,2)<=nfiles_use(iset)-min_present);
    nglomeruli_use(iset)=length(glomeruli_use{iset});
    resps_gu=resps_raw{iset}(:,glomeruli_use{iset},:);
    %
    %fill in missing data by affine interpolation
    % 
    resps_gur=reshape(resps_gu,[nstims(iset)*nglomeruli_use(iset),nfiles_use(iset)]);
    if ~exist('afalwt_opts') afalwt_opts=struct;end
    resps_tofill=isnan(resps_gur);
    if_canfill=1;
    if any(all(resps_tofill==1,1)) % JDD 10/8 some file has no data.
        disp('cannot fill in missing data, no stimuli present for some (glomerulus,file) pair')
        if_canfill=0;
    end
    if any(all(resps_tofill==1,2)) % some stim/glom pair is never present.
        disp('cannot fill in missing data, no (glomerulus,file) pair present for some stimulus')
        if_canfill=0;
    end
    if if_canfill
        %error('At canfill')
        [afalwt_fit,afalwt_b_change,afalwt_optsused]=afalwt(resps_gur,1-resps_tofill,afalwt_opts);
        resps_gur_fitted=(afalwt_fit.x_true*afalwt_fit.b_norm+repmat(afalwt_fit.a,size(resps_gur,1),1)); %interpolated data
        resps_gur_filled=resps_gur;
        resps_gur_filled(resps_tofill)=resps_gur_fitted(resps_tofill);
        resps_gu_filled=reshape(resps_gur_filled,[nstims(iset) nglomeruli_use(iset) nfiles_use(iset)]);
    else
        %error('not canfill')
        resps_gu_filled=resps_gu;
    end
    %
    disp(sprintf('%4.0f NaN values filled in',sum(resps_tofill(:))));
    %
    %show data, before glomerulus selection and after selection, normalized within each dataset
    %
    for ifig=1:3
        switch ifig
            case 1
                figname='all raw data';
                resps_plot=resps_raw{iset};
                glomeruli_plot=[1:nglomeruli];
            case 2
                figname='raw data from selected glomeruli';
                resps_plot=resps_gu;
                glomeruli_plot=glomeruli_use{iset};
            case 3
                figname='raw data from selected glomeruli with missing data filled in';
                resps_plot=resps_gu_filled;
                glomeruli_plot=glomeruli_use{iset};
        end
        figname=cat(2,sprintf('set %2.0f: ',iset),figname);
        figure;
        set(gcf,'Position',[50 100 1800 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figname);
        [nr,nc]=nicesubp(nfiles_use(iset));
        for ifile_ptr=1:nfiles_use(iset)
            ifile=files_use{iset}(ifile_ptr);
            subplot(nr,nc,ifile_ptr);
            imagesc(resps_plot(:,:,ifile_ptr));
            minmax=[min(min(resps_plot(:,:,ifile_ptr),[],'omitnan')),max(max(resps_plot(:,:,ifile_ptr),[],'omitnan'))];
            title_string=sprintf('set %1.0f file %2.0f: %s  [%6.3f %6.3f]',iset,ifile,strrep(filenames_short{iset}{ifile},'.mat',''),minmax);
            title(title_string,'Interpreter','none');
            set(gca,'FontSize',7);
            set(gca,'XTick',[1:length(glomeruli_plot)]);
            set(gca,'XTickLabel',glomeruli(glomeruli_plot));
            set(gca,'YTick',[1:nstims(iset)]);
            set(gca,'YTickLabel',stim_labels_set{iset});
        end
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,figname,'Interpreter','none');
        axis off;
    end %ifig
    %
    if any(isnan(resps_gu_filled(:)))
        error('Cannot proceed. Not all NaNs have been filled in.')
    end
    stim_labels_set{iset}=stimulus_names_set{iset};
    resps_set{iset}=reshape(afalwt_fit.x_true,[nstims(iset) nglomeruli_use(iset)]); %use regression slope as response measure
    if if_restore_size %restorr size if needed
        resps_set{iset}=resps_set{iset}*geomean(afalwt_fit.b_norm);
    end
    %
    %condition the data and stimulus names
    %
    if (if_submean)
        resps_set{iset}=resps_set{iset}-repmat(mean(resps_set{iset},1),nstims(iset),1);
    end
    %
    for istim=1:nstims(iset)
        if contains(stim_labels_set{iset}{istim},'@')
            stim_labels_set{iset}{istim}=stim_labels_set{iset}{istim}(1:min(find(stim_labels_set{iset}{istim}=='@')-1));
        end
        stim_labels_set{iset}{istim}=deblank(stim_labels_set{iset}{istim});
    end
end
%
%calculations common to all datasets
%
resp_range=[Inf -Inf];
for iset=1:nsets
    disp(sprintf('set %2.0f has %3.0f files used (of %3.0f), %3.0f stimuli, and %3.0f glomeruli used',iset,nfiles_use(iset),nfiles(iset),nstims(iset),nglomeruli_use(iset)))
    resp_range(1)=min(resp_range(1),min(resps_set{iset}(:)));
    resp_range(2)=max(resp_range(2),max(resps_set{iset}(:)));
end
dsid_show=sprintf('affine merging %2.0f files: %s,...,%s',sum(nfiles),deblank(dsid{1}(files_use{1}(1),:)),deblank(dsid{nsets}(files_use{nsets}(nfiles_use(nsets)),:)));
%
%overall response histogram and quantiles
%
figure;
set(gcf,'Position',[50 100 800 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','histogram');
for iset=1:nsets
    disp(sprintf('quantiles for set %1.0f',iset))
    subplot(2,nsets,iset);
    hist(resps_set{iset}(:),hist_bins);
    xlabel('response')
    ylabel('counts');
    set(gca,'XLim',resp_range);
    title(sprintf('set %1.0f',iset))
    quantiles=quantile(resps_set{iset}(:),hist_quantiles);
    for k=1:length(hist_quantiles)
        axes('Position',[0.01+(iset-1)/nsets,0.06+0.04*k,0.01,0.01]); %for text
        qt=sprintf(' quantile %5.3f: %7.3f',hist_quantiles(k),quantiles(k));
        disp(qt)
        text(0,0,qt,'Interpreter','none');
        axis off;
    end
end
axes('Position',[0.01,0.02,0.01,0.01]);
text(0,0,sprintf('restore size=%2.0f, subtract mean=%2.0f',if_restore_size,if_submean));
axis off;
%
%quantities calculated from all sets
%
stimulus_names=cell(0);
stim_labels_in=cell(0);
resp_range=[Inf -Inf];
for iset=1:nsets
    stimulus_names=[stimulus_names;stimulus_names_set{iset}];
    stim_labels_in=[stim_labels_in;stim_labels_set{iset}];
    resp_range=[min(resp_range(1),min(resps_set{iset}(:))) max(resp_range(2),max(resps_set{iset}(:)))];
end
% 
% stimulus_names_unique is a list of unique stimulus names, but in order of first occurrence in sets
% stim_labels_unique are the corresponding stim_labels, but does not contain concentration info
% (this is so that if the sets have non-overlapping names, then they will
% appear in the merged set in the original "stacked" order, rather than a merged alphabetical order
%
stimulus_names_unique=cell(0);
stim_labels_unique=cell(0);
for iset=1:nsets
    for istim=1:length(stimulus_names_set{iset})
        imatch=strmatch(stimulus_names_set{iset}(istim),stimulus_names_unique,'exact');
        if isempty(imatch)
            stimulus_names_unique{end+1,1}=stimulus_names_set{iset}{istim};
            stim_labels_unique{end+1,1}=stim_labels_set{iset}{istim};
        end
    end
end
if length(unique(stimulus_names))~=sum(nstims) | length(stimulus_names_unique)~=sum(nstims)
    warning('Not all stimulus names are unique');
end
if length(unique(stim_labels_in))~=sum(nstims) | length(stim_labels_unique)~=sum(nstims)
    warning('Not all stim labels are unique');
end
%
% determine contributions based on stimulus_names, since this also contains concentration
%
contribs=zeros(length(stimulus_names_unique),nsets);
for iset=1:nsets
    for istim=1:length(stimulus_names_set{iset})
        imatch=strmatch(stimulus_names_set{iset}{istim},stimulus_names_unique,'exact');
        if length(imatch)==1
            contribs(imatch,iset)=istim;
        end
    end
end
%
%initialize the quantities to save
%
f_base=struct;
for iset=1:nsets
    f_base.metadata{iset}=s{iset}{files_use{iset}(1)}.meta; %original metadata from Hong Lab
    f_base.dsid{iset}=dsid{iset}(files_use{iset},:); %data set ID, with special chars turned into -
end
f_base.resps=resps_set; %original responses
if if_restore_size==0
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in'; %original field for responses from Hong Lab
else
    f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in and restoring to original size';
end
%
%combined set either has intersection of glomeruli in components, or glomeruli in union, with missing responses set to 0
%
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
s_diag_set=cell(nsets,ncombm);
u_full_set=cell(nsets,ncombm);
v_full_set=cell(nsets,ncombm);
s_full_set=cell(nsets,ncombm);
%
%glomeruli_combined via intersection {1} and union {2}
glomeruli_combined=cell(1,ncombm);
glomeruli_combined{1}=[1:nglomeruli];
glomeruli_combined{2}=[];
for iset=1:nsets
    glomeruli_combined{1}=intersect(glomeruli_combined{1},glomeruli_use{iset});
    glomeruli_combined{2}=union(glomeruli_combined{2},glomeruli_use{iset});
end
%
%create concatenated datasets (stacked stimuli, just for display) and concatenated datasets 
%with multiple examples of each response averaged
%
if_ok=0;
while (if_ok==0)
    if_restrict_stims=getinp('1 to restrict stimuli to a subset before computing pcs of combined sets','d',[0 1]);
    if if_restrict_stims
        disp('**********'); %code from psg_coord_pipe_proc
        disp('Enter a file for the purposes of defining a stimulus set.')
        [sets_tn,ds_tn,sas_tn]=psg_get_coordsets(opts_read,[],[],1);
        subset_typenames=sas_tn{1}.typenames;
        disp('Stimulus set retrieved.');
        disp(subset_typenames');
        disp('**********');
        stims_combined_keep=[];
        for istim=1:length(stim_labels_unique)
            if length(strmatch(stim_labels_unique{istim},subset_typenames,'exact'))==1
                stims_combined_keep(end+1)=istim;
            end
        end   
        disp(sprintf('%3.0f stimuli of the original %3.0f in the combined files will be kept in the merged output',...
            length(stims_combined_keep),length(stim_labels_unique)));
        if_ok=getinp('1 if ok','d',[0 1]);
    else
        stims_combined_keep=[1:length(stim_labels_unique)];
        if_ok=1;
    end
end %if_ok
f_base.stimulus_names_orig=strvcat(stimulus_names_unique); %original stimulus names
f_base.stim_labels_orig=strvcat(stim_labels_unique); %shortened names for plotting
f_base.stimulus_names=f_base.stimulus_names_orig(stims_combined_keep,:);
f_base.stim_labels=f_base.stim_labels_orig(stims_combined_keep,:);
%
for icombm=1:2
    stims_sofar=0;
    resps_concat{icombm}=zeros(sum(nstims),length(glomeruli_combined{icombm}));
    resps_combined{icombm}=zeros(length(stimulus_names_unique),length(glomeruli_combined{icombm}));
    switch icombm %separate logic for intersection and union
        case 1             
            comb_label='intersection';
            comb_label_short='inter';
            for iset=1:nsets
                g_target=[1:length(glomeruli_combined{icombm})];
                g_source=find(ismember(glomeruli_use{iset},glomeruli_combined{icombm}));
                resps_concat{icombm}(stims_sofar+[1:nstims(iset)],g_target)=resps_set{iset}(:,g_source);
                stims_sofar=stims_sofar+nstims(iset);
                contribs_inds=find(contribs(:,iset)>0);
                resps_combined{icombm}(contribs_inds,g_target)=resps_combined{icombm}(contribs_inds,g_target)+resps_set{iset}(contribs(contribs_inds,iset),g_source);
            end
        case 2
            comb_label='union, missing set to 0';
            comb_label_short='union';
            for iset=1:nsets
                g_target=find(ismember(glomeruli_combined{icombm},glomeruli_use{iset}));
                g_source=[1:length(glomeruli_use{iset})];
                resps_concat{icombm}(stims_sofar+[1:nstims(iset)],g_target)=resps_set{iset}(:,g_source);
                stims_sofar=stims_sofar+nstims(iset);
                contribs_inds=find(contribs(:,iset)>0);
                resps_combined{icombm}(contribs_inds,g_target)=resps_combined{icombm}(contribs_inds,g_target)+resps_set{iset}(contribs(contribs_inds,iset),g_source);
            end
    end
    resps_combined{icombm}=resps_combined{icombm}./repmat(sum(contribs>0,2),1,length(glomeruli_combined{icombm})); %average responses according to number of sets in which they occur
%    maxdim=min(length(stimulus_names_unique),length(glomeruli_combined{icombm}))-if_submean;
    maxdim=min(length(stims_combined_keep),length(glomeruli_combined{icombm}))-if_submean;
    disp(sprintf(' combination method %20s: %3.0f responses from %2.0f glomeruli, max dim %2.0f',...
        comb_label,sum(nstims),length(glomeruli_combined{icombm}),maxdim))
    maxdim_use=maxdim;
    maxdim_allowed=maxdim;
    %
    %create coords by SVD and add metadata and plot
    %
    disp(sprintf('datasets merged by %s',comb_label));
    [f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f_base,resps_combined{icombm}(stims_combined_keep,:),maxdim,maxdim_use,if_submean);
    %save combined variables
    f_combm{icombm}=f;
    s_diag_all_combm{icombm}=s_diag_all;
    u_full_combm{icombm}=u_full;
    v_full_combm{icombm}=v_full;
    s_full_combm{icombm}=s_full;
    coords_all_combm{icombm}=coords_all;
    %
    f.glomeruli_names=glomeruli; %names of glomeruli in each dataset
    f.glomeruli_subset_use=glomeruli_use; %glomeruli used in each subset
    f.glomeruli_combined=glomeruli_combined{icombm};
    f.roi_names=glomeruli(glomeruli_use{icombm}); %for compatibilty with hlid_majaxes
    %
    %create coords by SVD for each subset, keeping all stimuli
    %
    stims_sofar=0;
    for iset=1:nsets
        resps_combined_set=resps_concat{icombm}(stims_sofar+[1:nstims(iset)],:);
        maxdim_set=min(size(resps_combined_set))-if_submean;
        disp(sprintf('component set %1.0f for merge by %s',iset,comb_label))
        [fset,s_diag_set{iset,icombm},u_full_set{iset,icombm},v_full_set{iset,icombm},s_full_set{iset,icombm}]=...
            hlid_coords_svd(struct(),resps_combined_set,maxdim_set,maxdim_set,if_submean);
        stims_sofar=stims_sofar+nstims(iset);
    end
    %
    if getinp(sprintf('1 if ok to write a coordinate file for %s',comb_label),'d',[0 1])
        data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid',cat(2,'merged_',comb_label_short));
        data_fullname_write_def=strrep(data_fullname_write_def,'kc_soma_nls','orn_merged');
        data_fullname_write_def=strrep(data_fullname_write_def,'odor17',cat(2,'odor',zpad(sum(nstims),3)));
        data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
        save(data_fullname_write,'-struct','f');
        disp(sprintf('wrote %s',data_fullname_write));
    end
    resps=resps_combined{icombm}(stims_combined_keep,:);
    roi_names=glomeruli(glomeruli_combined{icombm});
    stim_labels=f_base.stim_labels_orig(stims_combined_keep,:); %label only the stimuli we keep
    hlid_coords_plot;
    axes('Position',[0.5,0.02,0.01,0.01]); %for text
    text(0,0,comb_label,'Interpreter','none');
    axis off;
    %
    %special plots for combining datasets
    %
    figname=sprintf('combined sets via %s',comb_label);
    figure;
    set(gcf,'Position',[50 100 1800 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',figname);
    %full array of merged responses
    subplot(1,2,1)
    imagesc(resps_concat{icombm},resp_range);
    axis equal;
    axis tight;
    hold on;
    for iset=1:nsets-1
        stimct=sum(nstims(1:iset));
        plot([0 length(glomeruli_combined{icombm})]+[0.5 0.5],stimct+[0.5 0.5],'k','LineWidth',2);
    end
    title_string=figname;
    set(gca,'FontSize',7);
    set(gca,'XTick',[1:length(glomeruli_combined{icombm})]);
    set(gca,'XTickLabel',glomeruli(glomeruli_combined{icombm}));
    set(gca,'YTick',[1:sum(nstims)]);
    set(gca,'YTickLabel',stim_labels_unique);
    title(title_string,'Interpreter','none','FontSize',10);
    colorbar;
    %
    %compare PCs in individual sets and full set
    %
    ncols=2*nsets;
    nrows=2;
    for iset=1:nsets
        subplot(nrows,ncols,ncols/2+iset);
        dots=v_full_combm{icombm}'*v_full_set{iset,icombm};
        imagesc(abs(dots),[0 1]);
        xlabel(sprintf('eivec, set %1.0f',iset));
        ylabel('combined eivec');
        axis equal;
        axis tight;
        title('|dots|');
    end
    %
    %compare PCs in individual sets with each other
    %
    nset_pairs=nsets*(nsets-1)/2;
    ncols=2*max(2,nset_pairs);
    ipair=0;
    for jset=2:nsets
        for iset=1:jset-1
            ipair=ipair+1;
            subplot(2,ncols,ncols+ncols/2+ipair);
            dots=v_full_set{jset,icombm}'*v_full_set{iset,icombm};
            imagesc(abs(dots),[0 1]);
            xlabel(sprintf('eivec, set %1.0f',iset));
            ylabel(sprintf('eivec, set %1.0f',jset));
            axis equal;
            axis tight;
            title('|dots|');
        end
    end
    %
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,figname_raw,'Interpreter','none');
    axis off;
    axes('Position',[0.5,0.02,0.01,0.01]); %for text
    text(0,0,title_string,'Interpreter','none');
    axis off;
end %icombm
