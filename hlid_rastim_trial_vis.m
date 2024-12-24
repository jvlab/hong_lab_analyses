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
if ~exist('sub_labels') sub_labels={'',' (mean sub)'}; end %can replace  by a subset to shorten analysis
if ~exist('preproc_labels') preproc_labels={'raw','normalized'}; end 
npreproc=length(preproc_labels);
%
if ~exist('opts_plot') opts_plot=struct; end
%
nsubs=length(sub_labels);
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
if ~exist('dmax_def') dmax_def=10;end
dmax=getinp('max dimension of representational space to create','d',[1 nstims],dmax_def);
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
for ifile=1:nfiles
    for ipreproc=1:npreproc
        rm=resps_mean{ifile};
        rt=resps_trial{ifile};
        switch preproc_labels{ipreproc}
            case 'raw'
            case 'normalized'
                rm_norm=sqrt(sum(rm.^2,2));
                rm_norm(rm_norm==0)=1;
                rm=rm./repmat(rm_norm,[1 nrois_avail(ifile)]);
                rt_norm=sqrt(sum(rt.^2,3));
                rt_norm(rt_norm==0)=1;
                rt=rt./repmat(rt_norm,[1 1 nrois_avail(ifile)]);
        end
        nonans=find(all(~isnan(rm),2));
        npcs=max(dmax,min(length(nonans),size(rm,2)));
        [u_nonan,s,v]=svd(rm(nonans,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
        s=diag(s(1:npcs,1:npcs));
        u=nan(size(rm,1),npcs);
        u(nonans,:)=u_nonan(:,1:npcs);
        v=v(:,1:npcs);
        coords_mean=u*diag(s);
        %project rt into trial space, where coords are u*s
        coords_trial_nonan=zeros(size(rm,1),npcs,nrepts);
        for irept=1:nrepts
            coords_trial_nonan(:,:,irept)=reshape(rt(nonans,irept,:),length(nonans),size(rm,2))*v;
        end
        coords_trial=nan(size(rm,1),npcs,nrepts);
        coords_trial(nonans,:,:)=coords_trial_nonan;
        %
        tstring=sprintf('file %s: %s',dsids{1},preproc_labels{ipreproc});
        for icl=1:length(coord_lists)
            coord_list=coord_lists{icl};
            if max(coord_list(:))<=npcs
                figure;
                set(gcf,'Position',[100 100 1200 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',cat(2,tstring,sprintf(', dmax: %2.0f',max(coord_list(:)))));
                [nr,nc]=nicesubp(size(coord_list,1),0.7);
                for ic=1:size(coord_list,1)
                    dimlist=coord_list(ic,:);
                    opt_plot.axis_handle=subplot(nr,nc,ic);
                    [opts_plot_used,opts_plot_trial_used]=...
                        hlid_rastim_trial_vis_plot(coords_mean,coords_trial,dimlist,colors,stimulus_names_display,opts_plot);
                    if size(coord_list,1)==1 %add a dashed line to origin on large plots
                        for istim=1:nstims
                            plot3([0,coords_mean(istim,dimlist(1))],[0,coords_mean(istim,dimlist(2))],[0,coords_mean(istim,dimlist(3))],'k:');
                        end
                    end
                    hlid_rastim_trial_vis_legend;
                end
                axes('Position',[0.01,0.02,0.01,0.01]); %for text
                text(0,0,tstring,'Interpreter','none');
                axis off;
            end %max dim ok
        end %icl
    end %each preprocess
end %each file
