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
% 
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, SVD, HLID_RASTIM2COORDS_DEMO, HLID_POOL_PCACONTRIB.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
dim_text='dim'; %leadin for fields of d
%
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/data/kc_soma_nls/2022-10-10__fly01__megamat0.mat'; end
npool=getinp('number of files to pool','d',[1 Inf]);
%
metadata=cell(npool,1);
dsid=[];
resps=[];
nresps_each=zeros(1,npool);
for ipool=1:npool
    if_ok=0;
    while (if_ok==0)
        HongLab_fn=getinp(sprintf('Hong Lab file name for file %2.0f to be pooled',ipool),'s',[],HongLab_fn);
        %
        da=load(HongLab_fn);
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
    end %ifok
end %ipool
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
%do svd/pca
%
if (if_submean)
    resps_use=resps_use-repmat(mean(resps_use,1),nstims-length(stims_nan),1);
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
stim_labels=da.response_amplitude_stim.stim; %start with cell array
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
f.metadata=metadata; %original metadata from Hong Lab
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
f.coord_opts.aux.nresps_each=nresps_each; %number of responses from each dataset
f.resps=resps;
if getinp('1 if ok to write a coordinate file','d',[0 1])
    data_fullname_write_def=strrep(hlid_opts.coord_data_fullname_write_def,'dsid','pooled');
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
