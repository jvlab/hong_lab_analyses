%hlid_tnt_label_model_3c_compare: special-purpose script to compare tntlabel,
%tnt transformed by simplemodels,andsets tnt3c
%
% plots distance, and distance normalized by rms distance of reference dataset (tnt3c) from centroid
%
%
if ~exist('filepath') filepath='data/kc_tnt/'; end
if ~exist('filebase') filebase='hlid_odor17_coords_'; end
ndsets=6; %number of datasets
ndmax=6; %number of dimensions fit
d=cell(ndsets,1); %raw data
c=cell(ndmax,1); %coords
shortnames=cell(ndsets,1);
ref=6; %datset 6 is reference
nstims=17;
nr=3;
nc=6;
for idim=1:ndmax
    c{idim}=zeros(nstims,idim,ndsets);
end
dists=zeros(nstims+1,ndmax,ndsets);
d{1}.name='TNT-label';
d{1}.file=cat(2,filepath,filebase,'tnt-3c-model_procrustes-noscale.mat');
d{2}.name='uniform';
d{2}.file=cat(2,filepath,filebase,'tnt-3c-model_procrustes-scale.mat');
d{3}.name='affine';
d{3}.file=cat(2,filepath,filebase,'tnt-3c-model_affine-nooffset.mat');
d{4}.name='projective';
d{4}.file=cat(2,filepath,filebase,'tnt-3c-model_projective.mat');
d{5}.name='pwaffine';
d{5}.file=cat(2,filepath,filebase,'tnt-3c-model_pwaffine.mat');
d{6}.name='control';
d{6}.file=cat(2,filepath,filebase,'TNT3c_consensus_noscale.mat');
for idset=1:ndsets
    shortnames{idset}=d{idset}.name(1:3);
    d{idset}.data=load(d{idset}.file);
    for idim=1:ndmax
        vname=sprintf('dim%1.0f',idim);
        c{idim}(:,:,idset)=d{idset}.data.(vname);
    end
end
for idim=1:ndmax
    dists(1:nstims,idim,:)=sqrt(sum((c{idim}-repmat(c{idim}(:,:,ref),[1 1 ndsets])).^2,2)); %istim, idim, dataset(model)
end
%compute rms across stimuli
dists(end,:,:)=sqrt(mean(dists(1:nstims,:,:).^2.1));
stim_labels=strvcat(d{1}.data.stim_labels,'RMS');
%
models_show={[1:3],[1:5]};
for if_norm=0:1
    switch if_norm
        case 0
            ylab='dist';
        case 1
            ylab='norm dist';
            norm_vec=zeros(1,idim);
            for idim=1:ndmax
                vecs_ref=c{idim}(:,:,ref);
                vecs_ref=vecs_ref-repmat(mean(vecs_ref,1),nstims,1); %centered
                norm_vec(1,idim)=sqrt(sum(vecs_ref(:).^2)/nstims);
            end
    end
    for ims=1:length(models_show)
        tstring=[];
        mshow=models_show{ims};
        for k=1:length(mshow)
           idset=mshow(k);
           tstring=cat(2,tstring,' ',d{idset}.name);
        end
        tstring=cat(2,tstring,'-> control');
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        %
        dists_show=dists(:,:,mshow);

        for istim=1:nstims+1
            subplot(nr,nc,istim);
            resid_dist=reshape(dists_show(istim,:,:),[ndmax,length(mshow)])';
            if if_norm
                resid_dist_plot=resid_dist./repmat(norm_vec,length(mshow),1);
                ymax=max(max(max(dists_show./repmat(norm_vec,[nstims+1,1,length(mshow)]))));
            else
                resid_dist_plot=resid_dist;
                ymax=max(dists_show(:));
            end
            plot(resid_dist_plot,'LineWidth',2); %plot each dimension
            set(gca,'XLim',0.5+[0 length(mshow)]);
            set(gca,'XTick',[1:length(mshow)]);
            set(gca,'XTickLabel',shortnames(mshow));
            ylabel(ylab);
            set(gca,'YLim',[0 ymax]);
            title(stim_labels(istim,:));
            if (istim==nstims+1)
                legend(num2str([1:ndmax]'),'Location','NorthEast');
            end
        end
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end
end
