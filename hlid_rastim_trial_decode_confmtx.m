
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
%
% if_eachsubsamp: 1 to plot each subsample, 0 to plot average across subsamples, defaults to 0
%
%  See also: HLID_RASTIM_TRIAL_DECODE.
%
nsubsamps_use=size(results.subsamps_list_used,1);
if ~exist('if_eachsubsamp') if_eachsubsamp=0; end
if if_eachsubsamp
    subsamp_range=[1 nsubsamps_use];
else
    subsamp_range=[0 0];
end
nrows=max(results.dmax-results.dmin+1,2);
ncols=max(results.ndecs,3);
for isub=1:results.nsubs
    for ipreproc=1:results.npreprocs
        for isubsamp=subsamp_range(1):subsamp_range(2)
            fcs=zeros(results.dmax,results.ndecs);
            %plot heatmaps
            figure;
            set(gcf,'Position',[100 50 1400 950]);
            for id=results.dmin:results.dmax
                irow=id-results.dmin+1;
                for idec=1:results.ndecs
                    icol=idec;
                    subplot(nrows,ncols,icol+(irow-1)*ncols);
                    % confusion_matrices=zeros(nstims,nstims,ndecs,dmax,nsubs,npreprocs,nsubsamps_use); % d1: actual stim, d2: decoded stim, d3: decision rule, d4: dmax, d5: sub mean d6: normalize, d7: subsample set
                    if isubsamp>0
                        confmtx=results.confusion_matrices(:,:,icol,id,isub,ipreproc,isubsamp);
                    else
                        confmtx=sum(results.confusion_matrices(:,:,icol,id,isub,ipreproc,:),7);
                    end
                    fcs(id,idec)=sum(diag(confmtx))/sum(confmtx(:));
                    imagesc(confmtx);
                    set(gca,'XTick',[1:results.nstims]);
                    set(gca,'XTickLabel',results.stimulus_names_display);
                    set(gca,'YTick',[1:results.nstims]);
                    set(gca,'YTickLabel',results.stimulus_names_display);
                    axis square;
                    colormap hot;
                    title(sprintf('dim%1.0f [%s] fc %5.3f',id,results.dec_labels{idec},fcs(id,idec)),'Interpreter','none');
                end
            end
            %
            axes('Position',[0.01,0.05,0.01,0.01]); %for text
            if if_eachsubsamp
                text(0,0,sprintf('subsamp %2.0f of %2.0f',isubsamp,nsubsamps_use),'Interpreter','none');
            else
                text(0,0,sprintf('all subsamps of %2.0f',nsubsamps_use),'Interpreter','none');
            end
            axis off;
            %
            axes('Position',[0.01,0.03,0.01,0.01]); %for text
            text(0,0,sprintf('%s %s all folds, %s to %s',results.sub_labels{isub},results.preproc_labels{ipreproc},results.dsids{1},results.dsids{end}),'Interpreter','none');
            axis off;
            %
            axes('Position',[0.01,0.01,0.01,0.01]); %for text
            text(0,0,results.xv_label,'Interpreter','none');
            axis off;
            %
            %plot fraction correct as line graphs
            %
            figure;
            set(gcf,'Position',[100 100 1200 800]);
            %
            plot([1:results.dmax],fcs,'LineWidth',2);
            xlabel('dimension');
            set(gca,'XTick',[1:results.dmax]);
            set(gca,'XTickLabel',[1:results.dmax]);
            set(gca,'YLim',[0 max(fcs(:))]);
            ylabel('fraction correct');
            legend(results.dec_labels,'Location','best');
            %
            axes('Position',[0.01,0.05,0.01,0.01]); %for text
            if if_eachsubsamp
                text(0,0,sprintf('subsamp %2.0f of %2.0f',isubsamp,nsubsamps_use),'Interpreter','none');
            else
                text(0,0,sprintf('all subsamps of %2.0f',nsubsamps_use),'Interpreter','none');
            end
            axis off;
            %
            axes('Position',[0.01,0.03,0.01,0.01]); %for text
            text(0,0,sprintf('%s %s all folds, %s to %s',results.sub_labels{isub},results.preproc_labels{ipreproc},results.dsids{1},results.dsids{end}),'Interpreter','none');
            axis off;
            %
            axes('Position',[0.01,0.01,0.01,0.01]); %for text
            text(0,0,results.xv_label,'Interpreter','none');
            axis off;

        end %isubsamp
    end %ipreproc
end %isub
