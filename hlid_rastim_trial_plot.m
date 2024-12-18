% hlid_rastim_trial_plot: utility script for hlid_rastim_trial_pca
% to plot heatmap of either stimulus-averaged responses or single trial responses
plot_pos=cell(1,2);
plot_ptr=cell(1,2);
mxy=ones(1,2);
dxy=ones(1,2);
ticks_xy=cell(1,2);
for ixy=1:2
    switch size(dr,ixy)
        case nstims %stimulus-averaged
            plot_pos{ixy}=plot_pos_stim;
            plot_ptr{ixy}=ptr_stim;
            mxy(ixy)=1;
            dxy(ixy)=nstims;
            ticks_xy{ixy}=ticks_stim;
        case nstims*nrepts %trial-averaged
            plot_pos{ixy}=plot_pos_trial;
            plot_ptr{ixy}=ptr_trial;
            mxy(ixy)=nrepts;
            dxy(ixy)=nstims*nrepts;
            ticks_xy{ixy}=ticks_trial;
    end
end
d=zeros(dxy);
d(plot_pos{1}(plot_ptr{1}),plot_pos{2}(plot_ptr{2}))=dr(plot_ptr{1},plot_ptr{2});
imagesc(d'); %transpose so that first dim is horizontal
hold on;
for k=1:nstims-1
    plot(0.5+mxy(1)*[k,k],0.5+mxy(2)*[0 nstims],'k');
    plot(0.5+mxy(1)*[0 nstims],0.5+mxy(2)*[k,k],'k');
end
axis square;
set(gca,'Xtick',ticks_xy{1});
set(gca,'YTick',ticks_xy{2});
set(gca,'XTickLabel',display_order);
set(gca,'YTickLabel',display_order);
