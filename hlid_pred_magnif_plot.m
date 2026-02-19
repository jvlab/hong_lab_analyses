%hlid_pred_magnif_plot: plot results of hlid_pred_magnif_demo
%
%   See also:  HLID_PRED_MAGNIF_DEMO,
%   RS_IMPORT_COORDSETS, RS_KNIT_COORDSETS, RS_GEOFIT, RS_XFORM_APPLY, COOTODSQ.
%
hlid_setup;
% results.filenames_orn=filenames_short_orn;
% results.filenames_kc=filenames_short_kc;
% results.nstims=nstims;
% results.stim_labels=stim_labels;
% results.if_submean=if_submean;
% results.if_restore_size=if_restore_size;
% results.files_use_orn=files_use_orn;
% results.files_use_kc=files_use_kc;
% results.min_present=min_present;
% results.dim_max=dim_max;
% results.drop_list=drop_list;
% results.drop_stim_max=drop_stim_max;
% results.dists=dists;
% results.dists_kc_out_of_sample=dists_kc_out_of_sample;
% results.transforms=transforms;
%
dimlist_plot=getinp('dimension list to plot','d',[1 results.dim_max],[2:min(7,results.dim_max)]);
ndims_plot=length(dimlist_plot);
%
%scattergrams of predicted distance vs actual distance
%
if ~exist('dist_range') dist_range=100; end
model_names={'procrustes','affine'};
nmodels=length(model_names);
utria_sel=find(triu(ones(results.nstims),1)==1); %select values in upper triangular part of a matrix of size nstims
for inout=1:2 %insample vs out of sample
    dists_models=cell(1,nmodels);
    switch inout
        case 1
            inout_label='in-sample';
            modelvals=results.dists{1};
            for im=1:nmodels
                dists_models{im}=modelvals.(cat(2,'kc_',model_names{im}));
            end
        case 2
            inout_label='out-of-sample';
            modelvals=results.dists_kc_out_of_sample;
            for im=1:nmodels
                dists_models{im}=modelvals.(model_names{im});
            end
    end
    dists_data=results.dists{1}.kc_data;
    tstring=sprintf('kc distances, data vs model, %s',inout_label);
    figure;
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'Numbertitle','off');
    set(gcf,'Name',tstring);
    %
    for im=1:nmodels
        for idim_ptr=1:ndims_plot
            idim=dimlist_plot(idim_ptr);
            subplot(nmodels,ndims_plot,idim_ptr+(im-1)*ndims_plot);
            xvals_all=dists_data(:,:,idim);
            yvals_all=dists_models{im}(:,:,idim);
            xvals=xvals_all(utria_sel);
            yvals=yvals_all(utria_sel);
            plot(xvals(:),yvals(:),'k.');hold on;
            plot([0 dist_range],[0 dist_range],'k:');
            %
            rmse=sqrt(mean(xvals(:)-yvals(:)).^2);
            text(0.1*dist_range,0.85*dist_range,sprintf('rmse=%5.1f',rmse));
            r=corr(xvals,yvals);
            text(0.1*dist_range,0.60*dist_range,sprintf('r=%5.3f',r));
            %
            title(sprintf('dim %2.0f',idim));
            xlabel('kc dist')
            ylabel(sprintf('pred: %s',model_names{im}));
            set(gca,'XLim',[0 dist_range]);
            set(gca,'YLim',[0 dist_range]);
            set(gca,'XTick',[0 0.5 1]*dist_range);
            set(gca,'YTick',[0 0.5 1]*dist_range);
            axis square;
        end
    end
    %
    axes('Position',[0.01,0.04,0.01,0.01]);
    text(0,0,cat(2,tstring,sprintf('; if submean %1.0f nstims %1.0f',results.if_submean,results.nstims)));
    axis off
    %
    axes('Position',[0.01,0.01,0.01,0.01]);
    text(0,0,sprintf('ORN: %s ... %s',results.filenames_orn{1},results.filenames_orn{end}),'Interpreter','none');
    axis off
    axes('Position',[0.50,0.01,0.01,0.01]);
    text(0,0,sprintf('KC: %s ... %s',results.filenames_kc{1},results.filenames_kc{end}),'Interpreter','none');
    axis off
end
