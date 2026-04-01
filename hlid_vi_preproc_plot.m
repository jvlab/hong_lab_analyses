% hlid_vi_preproc_plot: plot results from hlid_vi_preproc
%
% 20Mar26: convert from variance ratio to F ratio
% 31Mar26: add max_timepoints
%
%  See also:  HLID_VI_PREPROC.
%
if ~exist('eigs_show_list') eigs_show_list=[20 120]; end
if ~exist('logrange') logrange=10^2; end
if ~exist('n_fw_show') n_fw_show=n_fws; end
opts_read=filldefault(opts_read,'max_timepoints',0);
label_string=data_file;
if opts_read.max_timepoints>0
    label_string=cat(2,label_string,sprintf(' max timepoints: %2.0f',opts_read.max_timepoints));
end
if_done=0;
while(if_done==0)
    n_fw_show=getinp('number of filter levels to show (0 to end)','d',[0 n_fws],n_fw_show);
    if n_fw_show==0
        if_done=1;
    else
        for k=1:n_meas
            disp(sprintf(' %1.0f->%s',k,resp_measures{k}));
        end
        meas_show_list=getinp('measures to show','d',[1 n_meas],[1:n_meas]);
        n_meas_show=length(meas_show_list);
        eigs_show_list=getinp('eigenvalues to show','d',[1 n_repts*n_stims],eigs_show_list);
        logrange=getinp('log range','f',[1 Inf],logrange);
        %
        var_ratios_cum=cumsum(var_ratios_eacheiv_num.*eival_sqs)./cumsum(var_ratios_eacheiv_den.*eival_sqs); %variance ratios as eigenvalues are added in
        frat_factor=n_repts/(1-1/n_stims); %to convert hlid_varrats.ratio to frat
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
            for meas_show_ptr=1:n_meas_show
                meas_ptr=meas_show_list(meas_show_ptr);
                %scree plots
                subplot(n_meas_show,n_cols,1+(meas_show_ptr-1)*n_cols)
                semilogy(eival_sqs(1:n_eigs,1:n_fw_show,meas_ptr),'.-');
                xlabel('eigenvalue');
                ylabel(cat(2,resp_measures{meas_ptr},' var explained'));
                set(gca,'YLim',max(max(eival_sqs(:,:,meas_ptr)))*[1/logrange 1]);
                legend(fw_labels,'Location','best');
                %
                totvar=sum(eival_sqs(:,1:n_fw_show,meas_ptr),1);
                subplot(n_meas_show,n_cols,2+(meas_show_ptr-1)*n_cols)
                semilogy(eival_sqs(1:n_eigs,1:n_fw_show,meas_ptr)./repmat(totvar,n_eigs,1),'.-');
                xlabel('eigenvalue');
                ylabel(cat(2,resp_measures{meas_ptr},' frac var explained'));
                set(gca,'YLim',0.5*[1/logrange 1]);
                legend(fw_labels,'Location','best');
                %
                %variance ratio for each eiv
                %
                subplot(n_meas_show,n_cols,3+(meas_show_ptr-1)*n_cols);
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
                subplot(n_meas_show,n_cols,4+(meas_show_ptr-1)*n_cols);
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
                subplot(n_meas_show,n_cols,5+(meas_show_ptr-1)*n_cols);
                plot(fw_list,part_ratios(:,meas_ptr),'k.-');
                set(gca,'YLim',[1 n_repts*n_stims]);
                set(gca,'XTick',fw_list);
                xlabel('filtering width')
                ylabel('overall participation ratio');
                %
                %participation ratio for each eiv
                %
                subplot(n_meas_show,n_cols,6+(meas_show_ptr-1)*n_cols);
                plot(part_ratios_eacheiv(1:n_eigs,1:n_fw_show,meas_ptr),'.-');
                set(gca,'YLim',[1 min(n_repts,n_stims)]);
                set(gca,'YTick',[1:min(n_repts,n_stims)]);
                xlabel('eigenvalue')
                ylabel('indiv participation ratio (wts)');
            end %meas_ptr
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,label_string,'Interpreter','none');
            axis off
            %
            %composite scattergram of wt partipation ratios and F ratios
            %
            zlabels={'indiv part ratios (rep*stm)','indiv part ratios (spatem)'};
            zmax=[n_repts, resp_minlength];
            if exist('part_ratios_st')
                n_part=2;
            else
                n_part=1;
            end
            figure;
            set(gcf,'NumberTitle','off');
            set(gcf,'Position',[50 100 1100 800]);
            set(gcf,'Name',sprintf('eig characteristics: 1 to %1.0f',n_eigs));
            z=hlid_varrats(rand(1,n_repts,n_stims));
            for ipart=1:n_part %show either participation ratio for rept*stim weights, or for spatiotemp
                for meas_show_ptr=1:n_meas_show
                    meas_ptr=meas_show_list(meas_show_ptr);
                    subplot(n_meas_show,n_part,ipart+(meas_show_ptr-1)*n_part);
                    for fw_ptr=1:n_fw_show
                        if ipart==1
                            pr=part_ratios_eacheiv(1:n_eigs,fw_ptr,meas_ptr);
                        else
                            pr=part_ratios_st(1:n_eigs,fw_ptr,meas_ptr);
                        end
                        fr=frat_factor*var_ratios_eacheiv(1:n_eigs,fw_ptr,meas_ptr);
                        frp=1-fcdf(fr,z.fdof(1),z.fdof(2));
                        hs=plot3([1:n_eigs],frp,pr,'.','MarkerSize',10);
                        hold on;
                        set(hs,'Color',colors(mod(fw_ptr-1,size(colors,1))+1,:));
                    end
                    set(gca,'XLim',[0 n_eigs]);
                    set(gca,'YLim',[0 1]);
                    set(gca,'ZLim',[1 zmax(ipart)]); %largest possible range of participation ratio
                    if (ipart==1)
                        set(gca,'ZTick',[1:zmax(ipart)]);
                    end
                    grid on
                    box on
                    xlabel('eiv');
                    ylabel('p(F)');
                    zlabel(zlabels{ipart});
                    title(resp_measures{meas_ptr});
                    legend(fw_labels,'Location','Best');
                    set(gca,'View',[-15,10]);
                end
            end
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,label_string,'Interpreter','none');
            axis off
        end %eigs_show_list
     end %if_done
end
