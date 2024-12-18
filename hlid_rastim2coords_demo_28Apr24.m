%hlid_rastim2coords_demo: read calcium imaging data from Hong Lab, and convert 
% z-scored response amplitudes averaged within each stimulus, to a coordinate file suitable for the psg package
%
% z-scored response amplitudes are converted to coordinates by singular value decomposition.
% Optionally subtracting the mean is equivalent to principal components analysis.
% 
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, SVD.
%
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
dim_text='dim'; %leadin for fields of d
%
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/HongLabOrig_for_jdv/data/kc_soma_nls/2022-10-10__fly01__megamat0.mat'; end
%
HongLab_fn=getinp('Hong Lab file name','s',[],HongLab_fn);
%
da=load(HongLab_fn);
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
resps=da.response_amplitude_stim.mean_peak;
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
%do svd/pca
%
if (if_submean)
    resps_use=resps-repmat(mean(resps,1),nstims,1);
else
    resps_use=resps;
end
[u,s,v]=svd(resps_use); %resps_use=u*s*v', with u and v both orthogonal, so u*s=resps_use*v
s_diag_all=diag(s(1:maxdim_allowed,1:maxdim_allowed));
var_total=sum(s_diag_all);
s=s(1:maxdim,1:maxdim);
s_diag=diag(s);
u=u(:,1:maxdim);
v=v(:,1:maxdim);
disp('fraction of variance explained with each component')
disp(sprintf('%6.4f ',s_diag/var_total))
coords_all=u(:,1:maxdim)*s;
for idim=1:maxdim
    f.(cat(2,dim_text,sprintf('%1.0f',idim)))=coords_all(:,1:idim);
end
%
stim_labels=stimulus_names;
%
for istim=1:nstims
    if contains(stim_labels{istim},'@')
        stim_labels{istim}=stim_labels{istim}(1:min(find(stim_labels{istim}=='@')-1));
    end
    stim_labels{istim}=deblank(stim_labels{istim});
end
%
%create the fields to save
%
f.metadata=da.meta; %original metadata from Hong Lab
f.stimulus_names=strvcat(stimulus_names); %original stimulus names
f.dsid=dsid; %data set ID, with special chars turned into -
f.stim_labels=strvcat(stim_labels); %shortened names for plotting
f.coord_opts.method='svd';
f.coord_opts.resp_type='response_amplitude_stim'; %original field for responses from Hong Lab
f.coord_opts.maxdim=maxdim;
f.coord_opts.if_submean=if_submean;
f.coord_opts.aux.desc='resps, or resps-repmat(mean(resps,1),nstims,1), = u*s*transpose(v) ';
f.coord_opts.aux.u=u;
f.coord_opts.aux.s=s;
f.coord_opts.aux.v=v;
f.resps=resps;
if getinp('1 if ok to write a coordinate file','d',[0 1])
    data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid',dsid);
    data_fullname_write=getinp('coordinate file name to write','s',[],data_fullname_write_def);
    save(data_fullname_write,'-struct','f');
    disp(sprintf('wrote %s',data_fullname_write));
end
%
figname_raw=sprintf('raw responses: %s',dsid);
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Name',figname_raw);
set(gcf,'Position',[100 100 1200 800]);
%
subplot(3,1,1)
imagesc(resps);
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
plot([1:maxdim_allowed],s_diag_all/var_total,'k');
hold on;
plot([1:maxdim],s_diag_all(1:maxdim)/var_total,'k*');
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
