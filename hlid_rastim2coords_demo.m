%hlid_rastim2coords_demo: read calcium imaging data from Hong Lab, and convert 
% z-scored response amplitudes averaged within each stimulus, to a coordinate file suitable for the psg package
%
% z-scored response amplitudes are converted to coordinates by singular value decomposition.
% SVD with subtracting the mean is equivalent to principal components analysis.
%
% Number of ROIs without NaNs must be at least as large as number of stimuli.
% 21May24: add removal of columns of NaN's from responses
% 15Sep24: modularize svd and partial creation of metadata; modularize plotting (in script)
% 06May25: pass stims_nan and stims_nonan to hlid_coords_svd
% 07May25: include a call to hlid_da_stimselect, to select stimuli and responses
%
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, SVD, HLID_RASTIM2COORDS_POOL, HLID_CSV2COORDS_DEMO,
%  HLID_COORDS_SVD, HLID_COORDS_PLOT, HLID_DA_STIMSELECT.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
if ~exist('opts_dasel')
    opts_dasel=struct;
end
%
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/data/kc_soma_nls/2022-10-10__fly01__megamat0.mat'; end
%
HongLab_fn=getinp('Hong Lab file name','s',[],HongLab_fn);
%
da=load(HongLab_fn);
[da,optsused_dasel]=hlid_da_stimselect(da,opts_dasel);
%
stimulus_names=da.response_amplitude_stim.stim';
dsid=da.meta.title;
%'-'can be used within fields of file name
dsid=strrep(dsid,'/','-');
dsid=strrep(dsid,'\','-');
dsid=strrep(dsid,'_','-');
while contains(dsid,'--')
    dsid=strrep(dsid,'--','-');
end
%
dsid=getinp('data set id','s',[],dsid);
nstims=length(stimulus_names);
%
resps_raw=da.response_amplitude_stim.mean_peak;
nancols=find(all(isnan(resps_raw),1));
disp(sprintf('ROIs found: %3.0f, ROIs that are all NaNs: %3.0f',size(resps_raw,2),length(nancols)));
resps=resps_raw(:,setdiff([1:size(resps_raw,2)],nancols));
disp(sprintf('ROIs analyzed: %3.0f',size(resps,2)));
disp(sprintf('unique stimuli found: %4.0f',nstims));
if size(resps,1)~=nstims
    warning(sprintf('number of responses (%2.0f) is not equal to number of stimuli (%2.0f)',size(resps,1),nstims));
end
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
f.metadata=da.meta; %original metadata from Hong Lab
f.stimulus_names=strvcat(stimulus_names); %original stimulus names
f.dsid=dsid; %data set ID, with special chars turned into -
f.stim_labels=strvcat(stim_labels); %shortened names for plotting
f.resps=resps; %original responses
f.coord_opts.resp_type='response_amplitude_stim'; %original field for responses from Hong Lab
%
%create coords by SVD and add metadata
%
[f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f,resps_use,maxdim,maxdim_use,if_submean,stims_nonan,stims_nan);
%
if getinp('1 if ok to write a coordinate file','d',[0 1])
    data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid',dsid);
    data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
    save(data_fullname_write,'-struct','f');
    disp(sprintf('wrote %s',data_fullname_write));
end
%plot
dsid_show=dsid;
hlid_coords_plot;
