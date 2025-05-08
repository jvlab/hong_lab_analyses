%hlid_orn_merge: read raw orn data files, handle missing data, merge the
%orn reponses, and write a coordinate file
%
% In contrast to hlid_rastim2coords_demo, hlid_rastim2coords_pool:
%  * this keeps track of responses of each ORN.
%  * Attempts to fill in missing data by fitting dataset to an additive and
% multiplicative constant for each dataset applied to a value that depends
% arbitrarily on stimulus and glomerulus, and does SVD on that -- this allows tracking of
% ORN's (as opposed to merging by Procrustes rotations of responses from individual datasets,
% or concatenation of ORN responses (as in hlad_rastim2coords_pool).
%
% 07May25: include a call to hlid_da_stimselect, to select stimuli and responses
%
%   See also:  HLID_SETUP, HLID_RASTM2COORDS_DEMO, HLID_RASTIM2COORDS_POOL, AFALWT, HLID_SVD_COORDS, HLID_PLOT_COORDS,
%   HLID_DA_STIMSELECT.
%
hlid_setup;
if ~exist('opts_dasel')
    opts_dasel=struct;
end
[filenames_short,pathname]=uigetfile('*fly*.mat','Select raw ORN data files','Multiselect','on');
nfiles=length(filenames_short);
s=cell(nfiles,1);
files_use=[];
nancols=cell(nfiles,1);
nanrows=cell(nfiles,1);
dsid=[]
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
resps=zeros(nstims,nglomeruli,nfiles_use);
glomeruli_missing=zeros(nglomeruli,nfiles_use);
for ifile_ptr=1:nfiles_use
    ifile=files_use(ifile_ptr);
    glomeruli_missing(nancols{ifile},ifile_ptr)=1;
    resps(:,:,ifile_ptr)=s{ifile}.response_amplitude_stim.mean_peak;
end
for ipresent=0:nfiles_use
    disp(sprintf('number of glomeruli present in at least %2.0f datasets: %3.0f (missing in %2.0f or less)',...
        ipresent,sum(sum(glomeruli_missing,2)<=(nfiles_use-ipresent)),nfiles_use-ipresent));
end
min_present=getinp('minimum number of preps that a glomerulus must be present in','d',[0 nfiles_use],1);
glomeruli_use=find(sum(glomeruli_missing,2)<=nfiles_use-min_present);
nglomeruli_use=length(glomeruli_use);
resps_gu=resps(:,glomeruli_use,:);
%
%fill in missing data by affine interpolation
% 
resps_gur=reshape(resps_gu,[nstims*nglomeruli_use,nfiles_use]);
if ~exist('afalwt_opts') afalwt_opts=struct;end
resps_tofill=isnan(resps_gur);
if_canfill=1;
if any(all(resps_tofill==1,1))
    disp('cannot fill in missing data, no stimuli present for some (glomerulus,file) pair')
    if_canfill=0;
end
if any(all(resps_tofill==1,2))
    disp('cannot fill in missing data, no (glomerulus,file) pair present for some stimulus')
    if_canfill=0;
end
if if_canfill
    [afalwt_fit,afalwt_b_change,afalwt_optsused]=afalwt(resps_gur,1-resps_tofill,afalwt_opts);
    resps_gur_fitted=(afalwt_fit.x_true*afalwt_fit.b_norm+repmat(afalwt_fit.a,size(resps_gur,1),1)); %interpolated data
    resps_gur_filled=resps_gur;
    resps_gur_filled(resps_tofill)=resps_gur_fitted(resps_tofill);
    resps_gu_filled=reshape(resps_gur_filled,[nstims nglomeruli_use nfiles_use]);
else
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
            resps_plot=resps;
            glomeruli_plot=[1:nglomeruli];
        case 2
            figname='raw data from selected glomeruli';
            resps_plot=resps_gu;
            glomeruli_plot=glomeruli_use;
        case 3
            figname='raw data from selected glomeruli with missing data filled in';
            resps_plot=resps_gu_filled;
            glomeruli_plot=glomeruli_use;
    end
    figure;
    set(gcf,'Position',[50 100 1800 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',figname);
    [nr,nc]=nicesubp(nfiles_use);
    for ifile_ptr=1:nfiles_use
        ifile=files_use(ifile_ptr);
        subplot(nr,nc,ifile_ptr);
        imagesc(resps_plot(:,:,ifile_ptr));
        minmax=[min(min(resps_plot(:,:,ifile_ptr),[],'omitnan')),max(max(resps_plot(:,:,ifile_ptr),[],'omitnan'))];
        title_string=sprintf('%s  [%6.3f %6.3f]',strrep(filenames_short{ifile},'.mat',''),minmax);
        title(title_string,'Interpreter','none');
        set(gca,'FontSize',7);
        set(gca,'XTick',[1:length(glomeruli_plot)]);
        set(gca,'XTickLabel',glomeruli(glomeruli_plot));
        set(gca,'YTick',[1:nstims]);
        set(gca,'YTickLabel',stim_labels);
    end
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,figname,'Interpreter','none');
    axis off;
end %ifig
%
if any(isnan(resps_gu_filled(:)))
    disp('Cannot proceed. Not all NaNs have been filled in.')
else
    resps=reshape(afalwt_fit.x_true,[nstims nglomeruli_use]); %use regression slope as response measure
    if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
    maxdim_allowed=min(size(resps))-if_submean;
    maxdim=getinp('maximum number of dimensions for coordinate file','d',[1 maxdim_allowed],maxdim_allowed);
    %
    %condition the data and stimulus names
    %
    if (if_submean)
        resps=resps-repmat(mean(resps,1),nstims,1);
    end
    maxdim_use=min(maxdim_allowed,size(resps,1));
    %
    stim_labels=stimulus_names;
    for istim=1:nstims
        if contains(stim_labels{istim},'@')
            stim_labels{istim}=stim_labels{istim}(1:min(find(stim_labels{istim}=='@')-1));
        end
        stim_labels{istim}=deblank(stim_labels{istim});
    end
    %
    %initialize the quantities to save
    %
    f=struct;
    f.metadata=s{files_use(1)}.meta; %original metadata from Hong Lab
    f.stimulus_names=strvcat(stimulus_names); %original stimulus names
    f.dsid=dsid(files_use,:); %data set ID, with special chars turned into -
    f.stim_labels=strvcat(stim_labels); %shortened names for plotting
    f.resps=resps; %original responses
    f.coord_opts.resp_type='response_amplitude_stim, after affine filling in'; %original field for responses from Hong Lab
    f.roi_names=glomeruli_check(glomeruli_use);
    %
    %create coords by SVD and add metadata
    %
    [f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f,resps,maxdim,maxdim_use,if_submean);
    dsid_show=sprintf('affine merging %2.0f files: %s,...,%s',nfiles,deblank(dsid(files_use(1),:)),deblank(dsid(files_use(nfiles_use),:)));
    %
    if getinp('1 if ok to write a coordinate file','d',[0 1])
        data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid','merged');
        data_fullname_write_def=strrep(data_fullname_write_def,'kc_soma_nls','orn_terminals');
        data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
        save(data_fullname_write,'-struct','f');
        disp(sprintf('wrote %s',data_fullname_write));
    end
    roi_names=f.roi_names;
    hlid_coords_plot;
end
