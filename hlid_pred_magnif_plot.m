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
if ~exist('dist_range_label') dist_range_label=100; end
if ~exist('ratio_range') ratio_range=[8 128]; end
if ~exist('dist_diff_range') dist_diff_range=[-20 70]; end
dist_range=max(dist_range_label,max(max(max(results.dists{1}.kc_data(:,:,dimlist_plot)))));
model_names={'procrustes','affine'};
nmodels=length(model_names);
utria_sel=find(triu(ones(results.nstims),1)==1); %select values in upper triangular part of a matrix of size nstims
%
%computed from kc distance
rmse=zeros(nmodels,ndims_plot,2); %root-mean-squared error, d1: model, d2:dim, d3: in-sample or out-of-sample
fvu=zeros(nmodels,ndims_plot,2); %fraction of variance unexplained, d1: model, d2:dim, d3: in-sample or out-of-sample
correl=zeros(nmodels,ndims_plot,2); %correlation, d1: model, d2:dim, d3: in-sample or out-of-sample
%computed from log of ratio of kc to orn distance
rmse_logratio=zeros(nmodels,ndims_plot,2); %root-mean-squared error, d1: model, d2:dim, d3: in-sample or out-of-sample
fvu_logratio=zeros(nmodels,ndims_plot,2); %fraction of variance unexplained, d1: model, d2:dim, d3: in-sample or out-of-sample
correl_logratio=zeros(nmodels,ndims_plot,2); %correlation, d1: model, d2:dim, d3: in-sample or out-of-sample
%
inout_labels={'in-sample','out-of-sample'};
for plot_type=1:4
    switch plot_type
        case 1
            plot_label='kc distances, data vs. model';
            m_low=1;
        case 2
            plot_label='kc/orn distances, data vs. model';
            m_low=1;
        case 3
            plot_label='kc distance heatmap';
            m_low=0;
        case 4
            plot_label='kc -model distance heatmap';
            m_low=1;
    end
    for inout=1:2 %insample vs out of sample
        dists_models=cell(1,nmodels);
        inout_label=inout_labels{inout};
        switch inout
            case 1
                 modelvals=results.dists{1};
                for im=1:nmodels
                    dists_models{im}=modelvals.(cat(2,'kc_',model_names{im}));
                end
            case 2
                modelvals=results.dists_kc_out_of_sample;
                for im=1:nmodels
                    dists_models{im}=modelvals.(model_names{im});
                end
        end
        dists_data=results.dists{1}.kc_data;
        tstring=sprintf('%s, %s',plot_label,inout_label);
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'Numbertitle','off');
        set(gcf,'Name',tstring);
        %
        for im=m_low:nmodels
            for idim_ptr=1:ndims_plot
                idim=dimlist_plot(idim_ptr);
                subplot(nmodels+1-m_low,ndims_plot,idim_ptr+(im-m_low)*ndims_plot);
                switch plot_type
                    case 1 % kc/orn dist
                        xvals_all=dists_data(:,:,idim);
                        yvals_all=dists_models{im}(:,:,idim);
                        xvals=xvals_all(utria_sel);
                        yvals=yvals_all(utria_sel);
                        plot(xvals,yvals,'k.');
                        hold on;
                        plot([0 dist_range],[0 dist_range],'k:');
                        %
                        rmse(im,idim_ptr,inout)=sqrt(mean((xvals-yvals).^2));
                        text(0.05*dist_range,0.95*dist_range,sprintf('rmse=%5.1f',rmse(im,idim_ptr,inout)),'FontSize',7);
                        fvu(im,idim_ptr,inout)=sum((xvals-yvals).^2)/sum((xvals-mean(xvals)).^2);
                        text(0.05*dist_range,0.85*dist_range,sprintf('fvu=%5.3f',fvu(im,idim_ptr,inout)),'FontSize',7);
                        correl(im,idim_ptr,inout)=corr(xvals,yvals);
                        text(0.05*dist_range,0.75*dist_range,sprintf('corr=%5.3f',correl(im,idim_ptr,inout)),'FontSize',7);
                        %
                        title(sprintf('dim %2.0f',idim));
                        xlabel('kc dist')
                        ylabel(sprintf('pred: %s',model_names{im}));
                        set(gca,'XLim',[0 dist_range]);
                        set(gca,'YLim',[0 dist_range]);
                        set(gca,'XTick',[0 0.5 1]*dist_range_label);
                        set(gca,'YTick',[0 0.5 1]*dist_range_label);
                        axis square;
                    case 2 %kc/orn dist
                        xvals_all=dists_data(:,:,idim)./results.dists{1}.orn_drop(:,:,idim);
                        yvals_all=dists_models{im}(:,:,idim)./results.dists{1}.orn_drop(:,:,idim);
                        xvals=log2(xvals_all(utria_sel));
                        yvals=log2(yvals_all(utria_sel));
                        plot(xvals,yvals,'k.');
                        hold on;
                        plot(log2(ratio_range),log2(ratio_range),'k:');
                        %
                        rmse_logratio(im,idim_ptr,inout)=sqrt(mean((xvals-yvals).^2));
                        text(log2(ratio_range(1)),0.05*log2(ratio_range(1))+0.95*log2(ratio_range(2)),sprintf('rmse=%5.1f',rmse_logratio(im,idim_ptr,inout)),'FontSize',7);
                        fvu_logratio(im,idim_ptr,inout)=sum((xvals-yvals).^2)/sum((xvals-mean(xvals)).^2);
                        text(log2(ratio_range(1)),0.15*log2(ratio_range(1))+0.85*log2(ratio_range(2)),sprintf('fvu=%5.3f',fvu_logratio(im,idim_ptr,inout)),'FontSize',7);
                        correl_logratio(im,idim_ptr,inout)=corr(xvals,yvals);
                        text(log2(ratio_range(1)),0.25*log2(ratio_range(1))+0.75*log2(ratio_range(2)),sprintf('corr=%5.3f',correl_logratio(im,idim_ptr,inout)),'FontSize',7);
                        %
                        title(sprintf('dim %2.0f',idim));
                        xlabel('kc/orn dist')
                        ylabel(sprintf('pred: %s',model_names{im}));
                        set(gca,'XLim',log2(ratio_range));
                        set(gca,'YLim',log2(ratio_range));
                        set(gca,'XTick',log2([ratio_range(1) geomean(ratio_range) ratio_range(2)]));
                        set(gca,'YTick',log2([ratio_range(1) geomean(ratio_range) ratio_range(2)]));
                        set(gca,'XTickLabel',[ratio_range(1) geomean(ratio_range) ratio_range(2)]);
                        set(gca,'YTickLabel',[ratio_range(1) geomean(ratio_range) ratio_range(2)]);
                        axis square;
                    case {3,4} %kc distance heatmap
                        if plot_type==3
                                if im==0
                                    dvals=dists_data(:,:,idim);
                                    mstring='data';
                                else
                                    dvals=dists_models{im}(:,:,idim);
                                    mstring=model_names{im};
                                end
                                zscale=[0 dist_range];
                        else
                            dvals=dists_data(:,:,idim)-dists_models{im}(:,:,idim);
                            mstring=sprintf('data - %s',model_names{im});
                            zscale=dist_diff_range;
                        end
                        imagesc(dvals,zscale);
                        set(gca,'XTick',[1:results.nstims]);
                        set(gca,'XTickLabel',results.stim_labels,'FontSize',7);
                        set(gca,'YTick',[1:results.nstims]);
                        set(gca,'YTickLabel',results.stim_labels,'FontSize',7);
                        title(sprintf(' dim %2.0f %s',idim,mstring));
                        axis square;
                end %plot_type
            end %idim_ptr
        end %im
        if ismember(plot_type,[3 4])
            axes('Position',[0.5,0.04,0.01,0.01]);
            text(0,0,sprintf('range: %5.2f to %5.2f',zscale));
            axis off
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
    end %inout
end %plot_type
%
%summary plots
%
figure;
set(gcf,'Position',[120 120 1000 800]);
set(gcf,'Numbertitle','off');
set(gcf,'Name','summary');
for if_logratio=1:2
    switch if_logratio
        case 1
            v_rmse=rmse;
            v_fvu=fvu;
            v_correl=correl;
            log_label='';
        case 2
            v_rmse=rmse_logratio;
            v_fvu=fvu_logratio;
            v_correl=correl_logratio;
            log_label=' (log2)';
    end
    subplot(3,2,if_logratio)
    plot(dimlist_plot,[v_rmse(:,:,1);v_rmse(:,:,2)]);
    set(gca,'ColorOrder',[1 0 0;1 0 0;0 0 1;0 0 1]);
    set(gca,'LineStyleOrder',{':','-',':','-'});
    set(gca,'LineStyleCyclingMethod','withcolor');
    xlabel('dim');
    legend({'proc in-samp','aff in-samp','proc out-of-samp','aff out-of-samp'},'Location','Best');
    title(cat(2,'rmse',log_label));
    %
    subplot(3,2,2+if_logratio)
    plot(dimlist_plot,[v_fvu(:,:,1);v_fvu(:,:,2)]);
    set(gca,'ColorOrder',[1 0 0;1 0 0;0 0 1;0 0 1]);
    set(gca,'LineStyleOrder',{':','-',':','-'});
    set(gca,'LineStyleCyclingMethod','withcolor');
    xlabel('dim');
    legend({'proc in-samp','aff in-samp','proc out-of-samp','aff out-of-samp'},'Location','Best');
    title(cat(2,'frac var unex',log_label));
    %
    subplot(3,2,4+if_logratio)
    plot(dimlist_plot,[v_correl(:,:,1);v_correl(:,:,2)]);
    set(gca,'ColorOrder',[1 0 0;1 0 0;0 0 1;0 0 1]);
    set(gca,'LineStyleOrder',{':','-',':','-'});
    set(gca,'LineStyleCyclingMethod','withcolor');
    xlabel('dim');
    set(gca,'YLim',[-0.1 1]);
    legend({'proc in-samp','aff in-samp','proc out-of-samp','aff out-of-samp'},'Location','Best');
    title(cat(2,'correl',log_label));
end
axes('Position',[0.01,0.04,0.01,0.01]);
text(0,0,sprintf('if submean %1.0f nstims %1.0f',results.if_submean,results.nstims));
axis off
%
axes('Position',[0.01,0.01,0.01,0.01]);
text(0,0,sprintf('ORN: %s ... %s',results.filenames_orn{1},results.filenames_orn{end}),'Interpreter','none');
axis off
axes('Position',[0.50,0.01,0.01,0.01]);
text(0,0,sprintf('KC: %s ... %s',results.filenames_kc{1},results.filenames_kc{end}),'Interpreter','none');
axis off
