%hlid_csv2coords_demo: read modeling data from Hong Lab in csv format, and convert 
% to a coordinate file suitable for the psg package
% 
% data (stimulus) columns are assumed to have @ in their header
%
% modeling data are converted to coordinates by singular value decomposition.
% SVD with subtracting the mean is equivalent to principal components analysis.
%
% 23Aug24: select data columns based on header, rather than assume 2:end
% 27Aug24: alphabetize the stimuli to match standard order in data files
% 15Sep24: modularize svd and partial creation of metadata; modularize plotting (in script)
%
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, SVD, HLID_RASTIM2COORDS_DEMO.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/MBmodels/megamat17_hemibrain_model_responses.csv'; end
%
HongLab_fn=getinp('Hong Lab csv file name','s',[],HongLab_fn);
%
data_cells=readcell(HongLab_fn);
stim_cols=find(contains(data_cells(1,:),'@')); %stimulus columns have @ in their header
resps_all=cell2mat(data_cells(2:end,stim_cols))';
stimulus_names=data_cells(1,stim_cols)';
%
fs_last=max([0 find(HongLab_fn=='/'),find(HongLab_fn=='\')]);
dsid=HongLab_fn(fs_last+1:end);
dsid=strrep(dsid,'_','-');
while contains(dsid,'--')
    dsid=strrep(dsid,'--','-');
end
%
dsid=getinp('data set id','s',[],dsid);
nstims=length(stimulus_names);
%
disp(sprintf('ROIs available: %3.0f',size(resps_all,2)));
resps_use_range=getinp('range to use','d',[1 size(resps_all,2)],[1 size(resps_all,2)]);
resps_use=resps_all(:,[resps_use_range(1):resps_use_range(2)]);
disp(sprintf('unique stimuli found: %4.0f',nstims));
if size(resps_use,1)>size(resps_use,2)
    warning(sprintf('number of responses (%2.0f) is greater than number of neural data channels or rois (%4.0f)',size(resps_use,1),size(resps_use,2)));
end
%
if_submean=getinp('1 to subtract mean across responses before creating coordinates','d',[0 1],0);
maxdim_allowed=min(size(resps_use))-if_submean;
maxdim=getinp('maximum number of dimensions for coordinate file','d',[1 maxdim_allowed],maxdim_allowed);
%
stim_labels=stimulus_names;
%
if_alphabetize=getinp('1 to alphabetize stimulus names to match data file order','d',[0 1],1);
for istim=1:nstims
    if contains(stim_labels{istim},'@')
        stim_labels{istim}=stim_labels{istim}(1:min(find(stim_labels{istim}=='@')-1));
    end
    stim_labels{istim}=deblank(stim_labels{istim});
end
disp('before sorting')
disp('stim_labels    stimulus_names')
for istim=1:nstims
    disp(sprintf('%15s %15s',stim_labels{istim},stimulus_names{istim}));
end
%
if if_alphabetize
    [stim_labels,sort_inds]=sortrows(stim_labels);
    stimulus_names=stimulus_names(sort_inds);
    resps_use=resps_use(sort_inds,:);
end
disp('after optional sort')
disp('stim_labels    stimulus_names')
for istim=1:nstims
    disp(sprintf('%15s %15s',stim_labels{istim},stimulus_names{istim}));
end
%
%condition the data
%
if (if_submean)
    resps_use=resps_use-repmat(mean(resps_use,1),nstims,1);
end
maxdim_use=min(maxdim_allowed,size(resps_use,1));
%
%initialize the quantities to save
%
f.metadata=struct(); %no metadata
f.stimulus_names=strvcat(stimulus_names); %stimulus names
f.dsid=dsid; %data set ID, with special chars turned into -
f.stim_labels=strvcat(stim_labels); %shortened names for plotting
f.resps=resps_use;
f.coord_opts.resp_type='simulation'; %original field for responses from Hong Lab
f.coord_opts.if_alphabetize=if_alphabetize;
%
%create coords by SVD and add metadata
%
[f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f,resps_use,maxdim,maxdim_use,if_submean);
%
if getinp('1 if ok to write a coordinate file','d',[0 1])
    data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid',dsid);
    data_fullname_write_def=strrep(data_fullname_write_def,'kc_soma_nls','models');
    data_fullname_write_def=strrep(data_fullname_write_def,'.csv','');
    data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
    save(data_fullname_write,'-struct','f');
    disp(sprintf('wrote %s',data_fullname_write));
end
%
dsid_show=dsid;
resps=resps_use; %no censoring in responses, this is a model
hlid_coords_plot;
