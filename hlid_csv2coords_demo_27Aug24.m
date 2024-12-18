%hlid_rastim2coords_demo: read modeling data from Hong Lab in csv format, and convert 
% to a coordinate file suitable for the psg package
% 
% data (stimulus) columns are assumed to have @ in their header
%
% modeling data are converted to coordinates by singular value decomposition.
% SVD with subtracting the mean is equivalent to principal components analysis.
%
% 23Aug24: select data columns based on header, rather than assume 2:end
% 27Aug24: alphabetize the stimuli to match standard order in data files
%
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, SVD, HLID_RASTIM2COORDS_DEMO.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
dim_text='dim'; %leadin for fields of d
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
%do svd/pca
%
if (if_submean)
    resps_use=resps_use-repmat(mean(resps_use,1),nstims,1);
end
maxdim_use=min(maxdim_allowed,size(resps_use,1));
[u,s,v]=svd(resps_use); %resps_use=u*s*v', with u and v both orthogonal, so u*s=resps_use*v
s_diag_all=diag(s(1:maxdim_use,1:maxdim_use));
var_total=sum(s_diag_all);
s=s(1:maxdim_use,1:maxdim_use);
u=u(:,1:maxdim_use);
v=v(:,1:maxdim_use);
if maxdim_use<maxdim %add zero eigenveectors and eigenvalues, should only happen if there are NaN stimuli
    u(stims_nonan,:)=u;
    u(stims_nan,:)=NaN;
    s(:,maxdim_use+1:maxdim)=0; %higher eigenvals and eigenvecs are zero
    s(maxdim_use+1:maxdim,:)=0;
    u(:,maxdim_use+1:maxdim)=0;
    v(:,maxdim_use+1:maxdim)=0;
end
%
s_diag=diag(s);
disp('fraction of variance explained with each component')
disp(sprintf('%6.4f ',s_diag/var_total))
coords_all=u*s;
for idim=1:maxdim
    f.(cat(2,dim_text,sprintf('%1.0f',idim)))=coords_all(:,1:idim);
end
%
%create the fields to save
%
f.metadata=struct(); %no metadata
f.stimulus_names=strvcat(stimulus_names); %stimulus names
f.dsid=dsid; %data set ID, with special chars turned into -
f.stim_labels=strvcat(stim_labels); %shortened names for plotting
f.coord_opts.method='svd';
f.coord_opts.resp_type='simulation'; %original field for responses from Hong Lab
f.coord_opts.maxdim=maxdim;
f.coord_opts.if_submean=if_submean;
f.coord_opts.if_alphabetize=if_alphabetize;
f.coord_opts.aux.desc='resps, or resps-repmat(mean(resps,1),nstims,1), = u*s*transpose(v) ';
f.coord_opts.aux.u=u;
f.coord_opts.aux.s=s;
f.coord_opts.aux.v=v;
f.resps=resps_use;
if getinp('1 if ok to write a coordinate file','d',[0 1])
    data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid',dsid);
    data_fullname_write_def=strrep(data_fullname_write_def,'kc_soma_nls','models');
    data_fullname_write_def=strrep(data_fullname_write_def,'.csv','');
    data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
    save(data_fullname_write,'-struct','f');
    disp(sprintf('wrote %s',data_fullname_write));
end
%
figname_raw=sprintf('raw responses: %s',dsid);
if if_submean
    figname_raw=cat(2,figname_raw,' mean subtracted');
else
    figname_raw=cat(2,figname_raw,' mean not subtracted');
end
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Name',figname_raw);
set(gcf,'Position',[100 100 1200 800]);
%
subplot(3,1,1)
imagesc(resps_use);
xlabel('roi')
ylabel('stimulus labels');
set(gca,'YTick',[1:nstims]);
set(gca,'YTickLabel',stim_labels);
colorbar;
%
subplot(3,1,2)
imagesc(v');
xlabel('roi')
ylabel('component');
set(gca,'YTick',[1:maxdim]);
colorbar;
%
subplot(3,3,7)
plot([1:length(s_diag_all)],s_diag_all/var_total,'k');
hold on;
plot([1:length(s_diag_all)],s_diag_all/var_total,'k*');
set(gca,'XLim',[0 maxdim_allowed+if_submean+0.5]);
set(gca,'XTick',[1:maxdim_allowed+if_submean]);
set(gca,'YLim',[0 1]);
xlabel('component');
ylabel('variance fraction');
title('variance explained by each component');
%
subplot(3,3,9);
imagesc(coords_all,[-1 1]*max(abs(coords_all(:))));
set(gca,'XTick',[1:maxdim]);
set(gca,'YTick',[1:nstims]);
set(gca,'YTickLabel',stim_labels);
title('coordinates')
colorbar;
%
axes('Position',[0.01,0.02,0.01,0.01]); %for text
text(0,0,figname_raw,'Interpreter','none');
axis off;
%
