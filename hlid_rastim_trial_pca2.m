%hlid_rastim_trial_pca2: auxiliary plots at single-trial level
%  
% Run this after hlid_rastim_trial_pca, or after loading results structure.
% See hlid_rastim_trial_pca for documentation and details.
%
% Computes and plots averages of distance heatmaps within trials, resids, orth.,
%  and plots scattergrams between stimulus-mean distances and trial-by-trial
%  distances, resids, orth.
% Only plots averages across files. 
%
%  See also:  HLID_RASTIM_TRIAL_PCA, HLID_RASTIM_TRIAL_PLOT.
%
%
%plot distances as heatmaps
%
if_signed_resit_orth=results.if_signed_resid_orth;
nrepts=results.nrepts;
nstims=results.nstims;
nfiles=results.nfiles;
dists=results.dists;
ndists=results.ndists;
dist_labels=results.dist_labels;
nprocs=results.nprocs;
proc_labels=results.proc_labels;
pcas=results.pcas;
npcas=results.npcas;
pca_labels=results.pca_labels;
nsubs=results.nsubs;
sub_labels=results.sub_labels;
metadata=results.metadata;
dsids=results.dsids;
stims_avail=results.stims_avail;
rois_avail=results.rois_avail;
rois=results.rois;
nprocs=results.nprocs;
%
display_order=results.display_order;
plot_pos_stim=results.plot_pos_stim;
ptr_stim=results.ptr_stim;
ticks_stim=results.ticks_stim;
plot_pos_trial=results.plot_pos_trial;
ptr_trial=results.ptr_trial;
ticks_trial=results.ticks_trial;
%
[nrows_dists,ncols_dists]=nicesubp(nprocs,0.7);

for isub=1:nsubs
    for idist=1:ndists
        for ifile=nfiles+1:nfiles+1 %only show averages across files
            if ifile<=nfiles
                tstring=sprintf('file %s, dist: %s %s',dsids{ifile},dist_labels{idist},sub_labels{isub});
            else
                tstring=sprintf('average across %2.0f files, %s to %s, dist: %s %s',nfiles,dsids{1},dsids{nfiles},dist_labels{idist},sub_labels{isub});
            end
            figure;
            set(gcf,'Position',[100 100 1400 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',tstring);
            dr_keep=cell(1,nprocs);
            for iproc=1:nprocs
                subplot(nrows_dists,ncols_dists,iproc);
                if (ifile<=nfiles)
                    dr=dists{ifile,idist,isub,iproc};
                else
                    dra=zeros(size(dists{1,idist,isub,iproc}));
                    for k=1:nfiles
                        dra(:,:,k)=dists{k,idist,isub,iproc};
                    end
                    dr=mean(dra,3,'omitnan');
                    %average across trials if needed
                    if size(dr,1)==nstims*nrepts
                        dr=reshape(mean(reshape(dr,[nrepts nstims size(dr,2)]),1),[nstims size(dr,2)]);
                    end
                    if size(dr,2)==nstims*nrepts
                        dr=reshape(mean(reshape(dr,[nstims nrepts nstims]),2),[nstims nstims]);
                    end
                end
                hlid_rastim_trial_plot;
                title(cat(2,proc_labels{iproc},sub_labels{isub}));
                dr_keep{iproc}=dr;
            end
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,tstring,'Interpreter','none');
            axis off;
            disp(sprintf('plotted %s',tstring));
            %
            %now plot scattergrams
            %
            figure;
            set(gcf,'Position',[100 100 1400 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'scattergrams ',tstring));
            for iproc=1:nprocs
                subplot(nrows_dists,ncols_dists,iproc);
                hg=plot(dr_keep{1}(:),dr_keep{iproc}(:),'k.');
                hold on;
                hd=plot(diag(dr_keep{1}),diag(dr_keep{iproc}),'r.');
                max_xy=max([get(gca,'XLim') get(gca,'Ylim')]);
                plot([0 max_xy],[0 max_xy],'k--');
                xlabel(proc_labels{1});               
                ylabel(proc_labels{iproc});
                set(gca,'XLim',[0 max(get(gca,'XLim'))]);
                set(gca,'YLim',[0 max(get(gca,'YLim'))]);
                title(cat(2,'dists ',sub_labels{isub}));               
                legend([hg;hd],{'off-diag','diag'},'Location','Best');
            end
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,tstring,'Interpreter','none');
            axis off;
        end %ifile
    end %idist
end %isub
