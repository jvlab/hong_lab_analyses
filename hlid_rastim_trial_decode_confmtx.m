
%hlid_rastim_trial_decode_confmtx: heatmaps of confusion matrices
% and line graphs of fraction correect
%
%runs on results structure from hlid_rastim_trial_decode
%
% heatmaps:
%* each row of plots is a dimension of the space, each col of plots is a decoding method
%* within each plot: row: a stimulus to be decoded, each col: what it is decoded as
%* one page for each choice of isub (mean sub or not) and ipreproc (normalize or not)
%* stimulus order is as in original data files
%
% line graphs:
%* fraction correct as a function of dimensionality, one line for each decoding method
%* transmitted information (raw and Treves-Panzeri debiased)
%
% if_eachsubsamp: 1 to plot confusion matrices and summary plot for each subsample;
%                 0 to plot confusion matrix and summary plot average across subsamples (default)
%                -1 as in 0, but superimposes individual subsamples on summary plot
% dmin_confmtx_show: minimum dimension to show confusoin matrix, defaults to results.dmin
% dmax_confmtx_show: maximum dimension to show confusion matrix, defaults to results.dmax
%
% 20Feb25: allow for dimension list to be non-contiguous
% 02Mar25: add information calculated from confusion matrix; add fraction correct and info in results
% 04Mar25: add plotting of analyses without embedding, present if results.dimlist contains Inf
% 06Mar25: fix bug with labeling of fraction correct; row-normalize confusion matrices and plot on scale of [0 1]
%
%  See also: HLID_RASTIM_TRIAL_DECODEm TBLXINFO_COUNT, TBLXTPBI.
%
%fill in, and backward compatibility
results=filldefault(results,'if_singleprep',0);
results=filldefault(results,'nsubsamps_avail',nchoosek(results.nfiles,results.nsets));
results=filldefault(results,'if_all_subsamps',double(size(results.subsamps_list_used,1)==results.nsubsamps_avail));
if ~isfield(results,'dimlist') %legacy: confusion matrix d4 goes from1 to results.dmax, but not all entries used
    dimlist=[1:results.dmax];
    dimlist_avail=results.dmin:results.dmax;
    dimlist_avail_ptr=dimlist_avail;
else %confusion matrix d4 is only the dimensions used
    dimlist=results.dimlist;
    dimlist_avail=dimlist;
    dimlist_avail_ptr=[1:length(dimlist)];
end
%
nsubsamps_use=size(results.subsamps_list_used,1);
%
if ~exist('dmin_confmtx_show') dmin_confmtx_show=results.dmin; end
dmin_confmtx_show=max(dmin_confmtx_show,results.dmin);
if ~exist('dmax_confmtx_show') dmax_confmtx_show=results.dmax; end
dmax_confmtx_show=min(dmax_confmtx_show,results.dmax);
%
if ~exist('if_eachsubsamp') if_eachsubsamp=0; end
%
if if_eachsubsamp>0
    subsamp_range=[1 nsubsamps_use];
else
    subsamp_range=[0 0];
end
dimlist_show=intersect([dmin_confmtx_show:dmax_confmtx_show],dimlist_avail);
if (max(dimlist)==Inf) %if there data for decoding without embedding
    dimlist_show=[dimlist_show Inf];
end
nrows=max(length(dimlist_show),2);
ncols=max(results.ndecs,3);
fcs=zeros(length(dimlist),results.ndecs,results.nsubs,results.npreprocs,1+nsubsamps_use); %d1: dim, d2: decode method, d3: sub mean or not, d4: norm or not, d5: 1+isubsamp
info_cm_raw=zeros(length(dimlist),results.ndecs,results.nsubs,results.npreprocs,1+nsubsamps_use); %d1: dim, d2: decode method, d3: sub mean or not, d4: norm or not, d5: 1+isubsamp
info_cm_debiased=zeros(length(dimlist),results.ndecs,results.nsubs,results.npreprocs,1+nsubsamps_use); %d1: dim, d2: decode method, d3: sub mean or not, d4: norm or not, d5: 1+isubsamp
% subsamp_label=cell(nsubsamps_use+1,1);
for isubsamp=0:nsubsamps_use
    for isub=1:results.nsubs
        for ipreproc=1:results.npreprocs
            for id=dimlist_avail
                id_ptr=find(dimlist==id);
                for idec=1:results.ndecs
                    icol=idec;
                    if isubsamp>0
                        confmtx=results.confusion_matrices(:,:,icol,id_ptr,isub,ipreproc,isubsamp);
                    else
                        confmtx=sum(results.confusion_matrices(:,:,icol,id_ptr,isub,ipreproc,:),7);
                    end
                    fcs(id_ptr,idec,isub,ipreproc,isubsamp+1)=sum(diag(confmtx))/sum(confmtx(:));
                    h=tblxinfo_count(confmtx);
                    info_cm_raw(id_ptr,idec,isub,ipreproc,isubsamp+1)=h;
                    info_cm_debiased(id_ptr,idec,isub,ipreproc,isubsamp+1)=h+tblxtpbi(confmtx,0); %Treves-Panzeri debiaser, all rows and columns
                end %idec
            end %id
        end %ipreproc
    end %isub
    if if_eachsubsamp>0
        subsamp_label{isubsamp+1}=sprintf('subsamp %2.0f of %2.0f',isubsamp,nsubsamps_use);
    else
        subsamp_label{isubsamp+1}=sprintf('all subsamps of %2.0f',nsubsamps_use);
    end
    if results.if_all_subsamps
        subsamp_label{isubsamp+1}=cat(2,subsamp_label{isubsamp+1},' (exhaustive)');
    end
    if results.if_singleprep
        subsamp_label{isubsamp+1}=cat(2,subsamp_label{isubsamp+1},' single prep mode');
    end
end %isubsamp
for isubsamp=subsamp_range(1):subsamp_range(2)
    for isub=1:results.nsubs
        for ipreproc=1:results.npreprocs
            %plot heatmaps
            label_proc=sprintf('%s %s',results.sub_labels{isub},results.preproc_labels{ipreproc});
            figure;
            set(gcf,'Position',[100 50 1400 950]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'confmtx: ',subsamp_label{isubsamp+1},' ',label_proc));
            for id=dimlist_show
                irow=find(dimlist_show==id);
                id_ptr=find(dimlist==id);
                for idec=1:results.ndecs
                    icol=idec;
                    if isubsamp>0
                        confmtx=results.confusion_matrices(:,:,icol,id_ptr,isub,ipreproc,isubsamp);
                    else
                        confmtx=sum(results.confusion_matrices(:,:,icol,id_ptr,isub,ipreproc,:),7);
                    end
                    %row-normalize
                    confmtx_rowsum=sum(confmtx,2);
                    confmtx_rowsum(confmtx_rowsum==0)=1;
                    confmtx_norm=confmtx./repmat(confmtx_rowsum,[1 size(confmtx,2)]);
                    subplot(nrows,ncols,icol+(irow-1)*ncols);
                    % confusion_matrices=zeros(nstims,nstims,ndecs,dmax,nsubs,npreprocs,nsubsamps_use); % d1: actual stim, d2: decoded stim, d3: decision rule, d4: dmax, d5: sub mean d6: normalize, d7: subsample set
                    imagesc(confmtx_norm,[0 1]);
                    set(gca,'XTick',[1:results.nstims]);
                    set(gca,'XTickLabel',results.stimulus_names_display);
                    set(gca,'YTick',[1:results.nstims]);
                    set(gca,'YTickLabel',results.stimulus_names_display);
                    axis square;
                    colormap hot;
                    if id==Inf
                        dimtag='no embed';
                    else
                        dimtag=sprintf('dim %1.0f',id);
                    end
                    % title(sprintf('%s [%s] fc %5.3f',dimtag,results.dec_labels{idec},fcs(id_ptr,idec)),'Interpreter','none'); %fix 06Mar25
                    title(sprintf('%s [%s] fc %5.3f',dimtag,results.dec_labels{idec},fcs(id_ptr,idec,isub,ipreproc,isubsamp+1)),'Interpreter','none');
                end
            end
            %
            axes('Position',[0.01,0.05,0.01,0.01]); %for text
            text(0,0,subsamp_label{isubsamp+1},'Interpreter','none');
            axis off;
            %
            axes('Position',[0.01,0.03,0.01,0.01]); %for text
            text(0,0,sprintf('%s %s all folds, %s to %s',results.sub_labels{isub},results.preproc_labels{ipreproc},results.dsids{1},results.dsids{end}),'Interpreter','none');
            axis off;
            %
            axes('Position',[0.01,0.01,0.01,0.01]); %for text
            text(0,0,results.xv_label,'Interpreter','none');
            axis off;
        end %ipreproc
    end %isub
end %isubsamp
%
%plot fraction correct as line graphs
%
dimlist_avail_plot=dimlist_avail;
if dimlist_avail(end)==Inf
    dimlist_avail_plot(end)=dimlist_avail_plot(end-1)+1;
    xticks=[1:results.dmax+1];
    xticklabels=strvcat(num2str([1:results.dmax]'),'no embed');
else
    xticks=[1:results.dmax];
    xticklabels=[1:results.dmax];
end
for isubsamp=subsamp_range(1):subsamp_range(2)
    figure;
    set(gcf,'Position',[100 100 1200 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'summary: ',subsamp_label{isubsamp+1}));
    for isub=1:results.nsubs
        for ipreproc=1:results.npreprocs
            subplot(results.nsubs,results.npreprocs,isub+(ipreproc-1)*results.nsubs);
            plot(dimlist_avail_plot,fcs(dimlist_avail_ptr,:,isub,ipreproc,isubsamp+1),'LineWidth',2);
            hold on;
            if if_eachsubsamp==-1
                color_order=get(gca,'ColorOrder');
                set(gca,'ColorOrder',color_order(1:results.ndecs,:)); %so that superimposed plots (if_eachsubsamp=-1) have same colors
                for iss=1:nsubsamps_use
                    plot(dimlist_avail_plot,fcs(dimlist_avail_ptr,:,isub,ipreproc,iss+1),':','LineWidth',1); %bug fix 04Mar25
                end
            end
            xlabel('dimension');
            set(gca,'XTick',xticks);
            set(gca,'XTickLabel',xticklabels);
            set(gca,'YLim',[0 max(fcs(:))]);
            ylabel('fraction correct');
            legend(results.dec_labels,'Location','best');
            title(sprintf('%s %s',results.sub_labels{isub},results.preproc_labels{ipreproc}),'Interpreter','none');
            %
        end %ipreproc
    end %isub
    axes('Position',[0.01,0.05,0.01,0.01]); %for text
    text(0,0,subsamp_label{isubsamp+1},'Interpreter','none');
    axis off;
    %
    axes('Position',[0.01,0.03,0.01,0.01]); %for text
    text(0,0,sprintf('all folds, %s to %s',results.dsids{1},results.dsids{end}),'Interpreter','none');
    axis off;
    %
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,results.xv_label,'Interpreter','none');
    axis off;
end %isubsamp
%
%plot info, raw and debiased,
%
ymax=max(max(info_cm_raw(:)),max(info_cm_debiased(:)));
for if_deb=0:1
    switch if_deb
        case 0
            tstring='raw info';
            info_plot=info_cm_raw;
        case 1
            tstring='debiased info';
            info_plot=info_cm_debiased;
    end
    for isubsamp=subsamp_range(1):subsamp_range(2)
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,tstring,' ',subsamp_label{isubsamp+1}));
        for isub=1:results.nsubs
            for ipreproc=1:results.npreprocs
                subplot(results.nsubs,results.npreprocs,isub+(ipreproc-1)*results.nsubs);
                plot(dimlist_avail_plot,info_plot(dimlist_avail_ptr,:,isub,ipreproc,isubsamp+1),'LineWidth',2);
                hold on;
                if if_eachsubsamp==-1
                    color_order=get(gca,'ColorOrder');
                    set(gca,'ColorOrder',color_order(1:results.ndecs,:)); %so that superimposed plots (if_eachsubsamp=-1) have same colors
                    for iss=1:nsubsamps_use
                        plot(dimlist_avail_plot,info_plot(dimlist_avail_ptr,:,isub,ipreproc,iss+1),':','LineWidth',1); %bug fix 04Mar25
                    end
                end
                xlabel('dimension');
                set(gca,'XTick',xticks);
                set(gca,'XTickLabel',xticklabels);
                set(gca,'YLim',[0 ymax]);
                ylabel(tstring);
                legend(results.dec_labels,'Location','best');
                title(sprintf('%s %s',results.sub_labels{isub},results.preproc_labels{ipreproc}),'Interpreter','none');
                %
            end %ipreproc
        end %isub
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,subsamp_label{isubsamp+1},'Interpreter','none');
        axis off;
        %
        axes('Position',[0.01,0.03,0.01,0.01]); %for text
        text(0,0,sprintf('all folds, %s to %s',results.dsids{1},results.dsids{end}),'Interpreter','none');
        axis off;
        %
        axes('Position',[0.01,0.01,0.01,0.01]); %for text
        text(0,0,results.xv_label,'Interpreter','none');
        axis off;
    end %isubsamp
end %if+deb

%save key calculations
results.fcs=fcs;
results.fcs_info_dims='d1: dim, d2: decode method, d3: sub mean or not, d4: norm or not, d5: 1+isubsamp';
results.info_cm_raw=info_cm_raw;
results.info_cm_debiased=info_cm_debiased;

