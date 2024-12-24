%hlid_rastim_trial_pca: visualize single-trial level data in representational space
%
% uses raw resonses (mean not subtracted), and normalized mean responses 
% (setting Euclidean distance to 1, so that distances are monotonically related to 1-correlation)
%
% psg_plotcoords used to enable lower-level control
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO,
%  HLID_RASTIM_TRIAL_READ, HLID_RASTIM_TRIAL_PCA, PSG_PLOTCOORDS,
%  PSG_VISUALIZE, HLID_RASTIM_TRIAL_VIS_PLOT, HLID_RASTIM_TRIAL_VIS_LEGEND.
%
hlid_setup;  %invoke hlid_localopts; set up opts_read and opts_plot
%
if ~exist('preproc_labels') preproc_labels={'raw','normalized'}; end %can replace by a subset to shorten analysis
npreproc=length(preproc_labels);
if ~exist('sub_labels') sub_labels={'',' (mean sub)'}; end %can replace by a subset to shorten analysis
nsubs=length(sub_labels);
if ~exist('st_labels') st_labels={'stims','trials'}; end %can replace by a subset to shorten analysis
nsts=length(st_labels);
%
if ~exist('opts_plot') opts_plot=struct; end
%
if ~exist('nrepts') nrepts=3; end %number of repeats
%
hlid_rastim_trial_read;
%
% metadata=cell(nfiles,1);
% dsids=cell(nfiles,1);
% resps_mean=cell(nfiles,1); mean responses within stimuli(nstims,nrois_avail)
% trial_ptrs=cell(nfiles,1); array of (nstims,nrepts)
% resps_trial=cell(nfiles,1); trial-by-trial responses, stimuli unscrambled, (nstims,nrepts,nrois_avail)
% trial_sequence=cell(nfiles,1); stimulus sequence, stimuli as strings
% stims_avail=cell(nfiles,1); list of available stimuli in each file, beginning at 1
% rois_avail=cell(nfiles,1); list of roi numbers kept for analysis, beginning at 1
% rois=cell(nfiles,1):  original rois
% nrois_avail(ifile): number of rois available
%
%create response space based on mean responses, and then project the
%single-trial responses into them
%
for k=1:nstims
    stimulus_names_display{k}=stimulus_names(k,1:-1+min(find(stimulus_names(k,:)==' ')));
end
%
if ~exist('opts_plot')
    opts_plot=struct;
end
if ~exist('colors') %colors for each trial
    colors={'r','m','c'};
end
%
if ~exist('coord_lists') %combinations of dimensions to show
    coord_lists{1}=[1 2 3];
    coord_lists{2}=[[1 2 3];[1 2 4];[1 3 4];[2 3 4]];
    coord_lists{3}=[2 3 4];
end
dmax=-inf;
for icl=1:length(coord_lists)
    dmax=max(dmax,max(coord_lists{icl}(:)));
end
dmax=getinp('max dimension of representational space to create','d',[dmax nstims],dmax);
%
for ifile=1:nfiles
    for isub=1:nsubs
        rs=resps_mean{ifile};
        rt=resps_trial{ifile};
        switch sub_labels{isub}
            case ''
            case ' (mean sub)'
            rs_xm=mean(rs,1,'omitnan'); %global mean-
            rs=rs-repmat(rs_xm,[nstims 1]);
            rt=rt-repmat(reshape(rs_xm,[1 1 nrois_avail(ifile)]),[nstims nrepts 1]);
        end
        for ipreproc=1:npreproc
            switch preproc_labels{ipreproc}
                case 'raw'
                case 'normalized'
                    rs_norm=sqrt(sum(rs.^2,2));
                    rs_norm(rs_norm==0)=1;
                    rs=rs./repmat(rs_norm,[1 nrois_avail(ifile)]);
                    rt_norm=sqrt(sum(rt.^2,3));
                    rt_norm(rt_norm==0)=1;
                    rt=rt./repmat(rt_norm,[1 1 nrois_avail(ifile)]);
            end
            for ist=1:nsts
                switch st_labels{ist}
                    case 'stims'
                        rpca=rs;
                    case 'trials'
                        rt_reshape=reshape(permute(rt,[2 1 3]),[nstims*nrepts,nrois_avail(ifile)]); %dim 1 is now stim (slow) and rept(fast)
                        rpca=rt_reshape;
                end
            nonans_pca=find(all(~isnan(rpca),2));
            npcs=min(dmax,min(length(nonans_pca),size(rpca,2)));
            [u_nonan,s,v]=svd(rpca(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
            s=diag(s(1:npcs,1:npcs));
            u=nan(size(rpca,1),npcs);
            u(nonans_pca,:)=u_nonan(:,1:npcs);
            v=v(:,1:npcs);
            %
            %coords_stim=u*diag(s) but compute it a different way so we can use the trial pca as well
            nonans_stim=find(all(~isnan(rs),2));
            coords_stim=nan(size(rs,1),npcs);
            coords_stim(nonans_stim,:)=rs(nonans_stim,:)*v;
            %project rt into pca space, where coords are u*s
            coords_trial_nonan=zeros(length(nonans_stim),npcs,nrepts);
            for irept=1:nrepts
                coords_trial_nonan(:,:,irept)=reshape(rt(nonans_stim,irept,:),length(nonans_stim),size(rs,2))*v;
            end
            coords_trial=nan(size(rs,1),npcs,nrepts);
            coords_trial(nonans_stim,:,:)=coords_trial_nonan;
            %
            tstring=sprintf('file %s: %s %s, pca on %s',dsids{1},sub_labels{isub},preproc_labels{ipreproc},st_labels{ist});
            for icl=1:length(coord_lists)
                coord_list=coord_lists{icl};
                if max(coord_list(:))<=npcs
                    figure;
                    set(gcf,'Position',[100 100 1200 800]);
                    set(gcf,'NumberTitle','off');
                    set(gcf,'Name',cat(2,tstring,sprintf(', dmax: %2.0f',max(coord_list(:)))));
                    [nrows,ncols]=nicesubp(max(size(coord_list,1),2),0.7); %leave a space for second plot
                    for ic=1:size(coord_list,1)
                        dimlist=coord_list(ic,:);
                        opt_plot.axis_handle=subplot(nrows,ncols,ic);
                        [opts_plot_used,opts_plot_trial_used]=...
                            hlid_rastim_trial_vis_plot(coords_stim,coords_trial,dimlist,colors,stimulus_names_display,opts_plot);
                        if size(coord_list,1)==1 %add a dashed line to origin on large plots
                            for istim=1:nstims
                                plot3([0,coords_stim(istim,dimlist(1))],[0,coords_stim(istim,dimlist(2))],[0,coords_stim(istim,dimlist(3))],'k:');
                            end
                            hlid_rastim_trial_vis_legend;
                            opts_plot.axis_handle=subplot(nrows,ncols,2);
                            %replot with mean response moved to origin
                            hlid_rastim_trial_vis_plot(zeros(size(coords_stim)),coords_trial-repmat(coords_stim,[1 1 nrepts]),dimlist,colors,{' '},...
                                setfield(opts_plot,'if_edges',1));
                            title('centered at mean resp');
                            hlid_rastim_trial_vis_legend;
                        else
                            hlid_rastim_trial_vis_legend;
                        end
                    end
                    axes('Position',[0.01,0.02,0.01,0.01]); %for text
                    text(0,0,tstring,'Interpreter','none');
                    axis off;
                end %max dim ok
            end %icl
            end %ists
        end %ipreproc
    end %isub
end %each file
