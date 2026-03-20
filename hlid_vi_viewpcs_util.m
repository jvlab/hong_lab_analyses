function hlid_vi_viewpcs_util(xyz,spatem_binned,if_hv,if_unifscale,bin_ranges,n_cols,col_off)
% hlid_vi_viewpcs_util(xyz,spatem_binned,if_hv,if_unifscale,bin_ranges,n_cols)
% is a utility to plot a volume as a column of slices
%
% xyz: a list of pixel positions, size(xyz,2)=3
% spatim_binned: data to plot, has dim 1=size(xyz,1), dim 2=size(bin_ranges,2) (already binned)
% if_hv: 1 for horiz, 2 for vertical (data transposed)
% if_unifscale: 1 for uniform scale, 0 for not
% bin_ranges: row 1, start bin, row 2: end bin
% n_cols: total columns in plot
% col_off: offset for first column
%
% See also: HLID_VI_VIEWPCS, HLID_VI_SPATIALFILTER.
%
font_size=7;
planes=unique(xyz(:,3));
n_planes=length(planes);
n_tbins=size(bin_ranges,2);
spatem_binned_max=max(abs(spatem_binned(:)));
for iplane=1:n_planes
    %code borrowed from hlid_vi_spatialfilter
    pxls_sel=find(xyz(:,3)==planes(iplane));
    xy_use=xyz(pxls_sel,[1:2]);
    xy_min=min(xy_use,[],1);
    xy_max=max(xy_use,[],1);
    xy_mod=xy_use-xy_min+1;
    data_mask=sparse(xy_mod(:,1),xy_mod(:,2),1);
    data_mask_full=full(data_mask);
    for tbin=1:n_tbins
        f=sparse(xy_mod(:,1),xy_mod(:,2),spatem_binned(pxls_sel,tbin));
        fu=full(f);
        fu(data_mask_full==0)=NaN; %mask out the NaN's
        subplot(n_planes,n_cols,col_off+tbin+(iplane-1)*n_cols);
        switch if_hv
            case 1
                if if_unifscale
                    imagesc(fu,spatem_binned_max*[-1 1]);
                else
                    imagesc(fu,max(abs(fu(:)))*[-1 1]);
                end
                xlabel(sprintf('frames [%2.0f %2.0f]',bin_ranges(:,tbin)),'FontSize',font_size);
                ylabel(sprintf('plane %2.0f',planes(iplane)),'FontSize',font_size);
                axis equal;
                set(gca,'YTick',1+[0 xy_max(1)-xy_min(1)]);
                set(gca,'YTickLabels',[xy_min(1),xy_max(1)],'FontSize',font_size);
                set(gca,'XTick',1+[0 xy_max(2)-xy_min(2)]);
                set(gca,'XTickLabels',[xy_min(2),xy_max(2)],'FontSize',font_size);
                axis tight;
            case 2
                if if_unifscale
                    imagesc(fu',spatem_binned_max*[-1 1]);
                else
                    imagesc(fu',max(abs(fu(:)))*[-1 1]);
                end
                xlabel(sprintf('frames [%2.0f %2.0f]',bin_ranges(:,tbin)),'FontSize',font_size);
                ylabel(sprintf('plane %2.0f',planes(iplane)),'FontSize',font_size);
                axis equal;
                set(gca,'XTick',1+[0 xy_max(1)-xy_min(1)]);
                set(gca,'XTickLabels',[xy_min(1),xy_max(1)],'FontSize',font_size);
                set(gca,'YTick',1+[0 xy_max(2)-xy_min(2)]);
                set(gca,'YTickLabels',[xy_min(2),xy_max(2)],'FontSize',font_size);
                axis tight;
        end %if_hv
    end %tbin
end
return
end
