% hlid_vi_preproc_plot: plot results from hlid_vi_preproc
%
% 20Mar26: convert from variance ratio to F ratio
%
%  See also:  HLID_VI_PREPROC.
%
if ~exist('eigs_show_list') eigs_show_list=[20 120]; end
if ~exist('logrange') logrange=10^2; end
if ~exist('n_fw_show') n_fw_show=n_fws; end
if_done=0;
while(if_done==0)
    n_fw_show=getinp('number of filter levels to show (0 to end)','d',[0 n_fws],n_fw_show);
    if n_fw_show==0
        if_done=1;
    else
        eigs_show_list=getinp('eigenvalues to show','d',[1 n_repts*n_stims],eigs_show_list);
        logrange=getinp('log range','f',[1 Inf],logrange);
        %
        var_ratios_cum=cumsum(var_ratios_eacheiv_num.*eival_sqs)./cumsum(var_ratios_eacheiv_den.*eival_sqs); %variance ratios as eigenvalues are added in
        frat_factor=n_repts/(1-1/n_stims); %to convert hlid_varrats.ratio to frt
        var_frats_cum=frat_factor*var_ratios_cum;
        %
        for eigs_show_ptr=1:length(eigs_show_list)
            n_eigs=min(eigs_show_list(eigs_show_ptr),n_repts*n_stims);
            fw_labels=cell(1,n_fw_show);
            for fw_ptr=1:n_fw_show
                fw_labels{fw_ptr}=sprintf('fw %2.0f',fw_list(fw_ptr));
            end
            figure;
            set(gcf,'NumberTitle','off');
            set(gcf,'Position',[50 50 1400 800]);
            set(gcf,'Name',sprintf('eigenvalue anlaysis: 1 to %1.0f',n_eigs));
            n_cols=6;
            for meas_ptr=1:n_meas
                %scree plots
                subplot(n_meas,n_cols,1+(meas_ptr-1)*n_cols)
                semilogy(eival_sqs(1:n_eigs,1:n_fw_show,meas_ptr),'.-');
                xlabel('eigenvalue');
                ylabel(cat(2,resp_measures{meas_ptr},' var explained'));
                set(gca,'YLim',max(max(eival_sqs(:,:,meas_ptr)))*[1/logrange 1]);
                legend(fw_labels,'Location','best');
                %
                totvar=sum(eival_sqs(:,1:n_fw_show,meas_ptr),1);
                subplot(n_meas,n_cols,2+(meas_ptr-1)*n_cols)
                semilogy(eival_sqs(1:n_eigs,1:n_fw_show,meas_ptr)./repmat(totvar,n_eigs,1),'.-');
                xlabel('eigenvalue');
                ylabel(cat(2,resp_measures{meas_ptr},' frac var explained'));
                set(gca,'YLim',0.5*[1/logrange 1]);
                legend(fw_labels,'Location','best');
                %
                %variance ratio for each eiv
                %
                subplot(n_meas,n_cols,3+(meas_ptr-1)*n_cols);
                semilogy(frat_factor*var_ratios_eacheiv(1:n_eigs,1:n_fw_show,meas_ptr),'.-');
                xlabel('eigenvalue')
                ylabel('indiv F ratio');
                hold on;
                %variance ratio for full data
                colors=get(gca,'ColorOrder');
                for fw_ptr=1:n_fw_show
                    hp=plot([1 n_eigs],repmat(var_frats(fw_ptr,meas_ptr),2,1));
                    set(hp,'Color',colors(mod(fw_ptr-1,size(colors,1))+1,:));
                end
                plot([1 n_eigs],[1 1],'k');
                ylim_vr=get(gca,'YLim');
                %
                %variance ratio, cumulative across eivs
                %
                subplot(n_meas,n_cols,4+(meas_ptr-1)*n_cols);
                semilogy(var_frats_cum(1:n_eigs,1:n_fw_show,meas_ptr),'.-');
                xlabel('eigenvalue')
                ylabel('cumul F ratio');
                hold on;
                %variance ratio for full data
                colors=get(gca,'ColorOrder');
                for fw_ptr=1:n_fw_show
                    hp=plot([1 n_eigs],repmat(var_frats(fw_ptr,meas_ptr),2,1));
                    set(hp,'Color',colors(mod(fw_ptr-1,size(colors,1))+1,:));
                end
                plot([1 n_eigs],[1 1],'k');
                set(gca,'YLim',ylim_vr);
                %
                %participation ratio for as function of filtering
                %
                subplot(n_meas,n_cols,5+(meas_ptr-1)*n_cols);
                plot(fw_list,part_ratios(:,meas_ptr),'k.-');
                set(gca,'YLim',[1 n_repts*n_stims]);
                set(gca,'XTick',fw_list);
                xlabel('filtering width')
                ylabel('overall participation ratio');
                %
                %participation ratio for each eiv
                %
                subplot(n_meas,n_cols,6+(meas_ptr-1)*n_cols);
                plot(part_ratios_eacheiv(1:n_eigs,1:n_fw_show,meas_ptr),'.-');
                set(gca,'YLim',[1 min(n_repts,n_stims)]);
                xlabel('eigenvalue')
                ylabel('indiv participation ratio');
            end
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,data_file,'Interpreter','none');
            axis off
        end %eigs_show_list
    end
end
