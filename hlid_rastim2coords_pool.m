%hlid_rastim2coords_pool: read calcium imaging data from Hong Lab, and convert 
% z-scored response amplitudes averaged within each stimulus, to a coordinate file suitable for the psg package
%
%  Differs from hlid_rastim2coords_demo in that several raw data files are pooled prior to SVD,
%  output file metadata is now a cell array, and dsid is a vertical concatenation of strings
% 
% z-scored response amplitudes are converted to coordinates by singular value decomposition.
% SVD with subtracting the mean is equivalent to principal components analysis.
%
% Number of ROIs without NaNs must be at least as large as number of stimuli.
%
% 02Jun24: added nresps_each, number of responses in each dataset, to coord_opts.aux
% 15Sep24: modularize svd and partial creation of metadata; modularize plotting (in script); fix plot label
% 30Oct24: use dialog box
% 06May25: pass stims_nan and stims_nonan to hlid_coords_svd
% 07May25: include a call to hlid_da_stimselect, to select stimuli and responses
%
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, SVD, HLID_RASTIM2COORDS_DEMO, HLID_POOL_PCACONTRIB, 
%  HLID_COORDS_SVD, HLID_COORDS_PLOT, HLID_DA_STIMSELECT.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/data/kc_soma_nls/2022-10-10__fly01__megamat0.mat'; end
if ~exist('opts_dasel')
    opts_dasel=struct;
end
npool_signed=0;
while npool_signed==0
    npool_signed=getinp('number of files to pool (<0 for dialog box)','d',[-Inf Inf]);
end
npool=abs(npool_signed);
%
metadata=cell(npool,1);
dsid=[];
resps=[];
nresps_each=zeros(1,npool);
if_ok=0;
ui_prompt='Select Hong lab files to be pooled';
ui_filter='_fly';
while (if_ok==0)
    if (npool_signed<0)
        [filenames_short,pathname,filter_index]=uigetfile(ui_filter,ui_prompt,'Multiselect','on');
    end
    %non-UI method
    for ipool=1:npool
        if (npool_signed>0)
            HongLab_fn=getinp(sprintf('Hong Lab file name for file %2.0f to be pooled',ipool),'s',[],HongLab_fn);
        else
            if ~iscell(filenames_short)
                filenames_short=cellstr(filenames_short);
            end
            disp(sprintf(' file %2.0f: %s',ipool,filenames_short{ipool}));
            HongLab_fn=cat(2,pathname,filesep,filenames_short{ipool});
        end
        da=load(HongLab_fn);
        [da,optsused_dasel]=hlid_da_stimselect(da,opts_dasel);
        metadata{ipool}=da.meta;
        stimulus_names_this=strvcat(da.response_amplitude_stim.stim');
        dsid_this=da.meta.title;
        %'-'can be used within fields of file name
        dsid_this=strrep(dsid_this,'/','-');
        dsid_this=strrep(dsid_this,'\','-');
        dsid_this=strrep(dsid_this,'_','-');
        while contains(dsid_this,'--')
            dsid_this=strrep(dsid_this,'--','-');
        end
        %
        dsid_this=getinp('data set id','s',[],dsid_this);
        dsid=strvcat(dsid,dsid_this);
        nstims_this=length(stimulus_names_this);
        if ipool==1
            nstims=nstims_this;
            stimulus_names=stimulus_names_this;
        end
        if_ok=1;
        if nstims~=nstims_this
            if_ok=0;
        else
            if any(stimulus_names~=stimulus_names_this)
                if_ok=0;
            end
        end
        %
        if (if_ok==0)
            disp('incompatible number of stimuli or stimulus names');
        else %file has compatible stimuli, so concatenate
            resps_raw=da.response_amplitude_stim.mean_peak;
            if size(resps_raw,1)~=nstims
                warning(sprintf('number of responses (%2.0f) is not equal to number of stimuli (%2.0f)',size(resps,1),nstims));
            end
            %
            nancols=find(all(isnan(resps_raw),1));
            disp(sprintf('ROIs found: %3.0f, ROIs that are all NaNs: %3.0f',size(resps_raw,2),length(nancols)));
            disp(sprintf('ROIs analyzed: %3.0f',size(resps_raw,2)));
            disp(sprintf('unique stimuli found: %4.0f',nstims));
            resps=[resps,resps_raw(:,setdiff([1:size(resps_raw,2)],nancols))];
            disp(sprintf('responses accumulated: %5.0f',size(resps,2)));
            nresps_each(ipool)=size(resps_raw,2)-length(nancols);
        end
    end %ipool
    if (if_ok)
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end %ifok
%
if size(resps,1)>size(resps,2)
    warning(sprintf('number of responses (%2.0f) is greater than number of neural data channels or rois (%4.0f)',size(resps,1),size(resps,2)));
end
%
if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
maxdim_allowed=min(size(resps))-if_submean;
maxdim=getinp('maximum number of dimensions for coordinate file','d',[1 maxdim_allowed],maxdim_allowed);
%
%find and temporarily remove stimuli that have nan responses.  if maxdim
%exceeds the number of remaining responses, the higher PCs and coords are zeros
% 
stims_nonan=find(all(~isnan(resps),2));
disp(sprintf('number of stimuli without NaN responses: %3.0f',length(stims_nonan)));
stims_nan=setdiff(1:nstims,stims_nonan);
resps_use=resps(stims_nonan,:);
%
%condition the data and stimulus names
%
if (if_submean)
    resps_use=resps_use-repmat(mean(resps_use,1),nstims-length(stims_nan),1);
end
maxdim_use=min(maxdim_allowed,size(resps_use,1));
%
stim_labels=da.response_amplitude_stim.stim; %start with cell array
for istim=1:nstims
    if contains(stim_labels{istim},'@')
        stim_labels{istim}=stim_labels{istim}(1:min(find(stim_labels{istim}=='@')-1));
    end
    stim_labels{istim}=deblank(stim_labels{istim});
end
%
%initialize the fields to save
%
f.metadata=metadata; %original metadata from Hong Lab
f.stimulus_names=strvcat(stimulus_names); %original stimulus names
f.dsid=dsid; %data set ID, with special chars turned into -
f.stim_labels=strvcat(stim_labels); %shortened names for plotting
f.resps=resps; %original responses
f.coord_opts.resp_type='response_amplitude_stim'; %original field for responses from Hong Lab
%
%create coords by SVD and add metadata
%
[f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f,resps_use,maxdim,maxdim_use,if_submean,stims_nonan,stims_nan);
f.coord_opts.aux.nresps_each=nresps_each; %number of responses from each dataset
%
if getinp('1 if ok to write a coordinate file','d',[0 1])
    data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid','pooled');
    data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
    save(data_fullname_write,'-struct','f');
    disp(sprintf('wrote %s',data_fullname_write));
end
%plot
dsid_show=sprintf('%2.0f files: %s,...,%s',npool,deblank(dsid(1,:)),deblank(dsid(end,:)));
hlid_coords_plot;
