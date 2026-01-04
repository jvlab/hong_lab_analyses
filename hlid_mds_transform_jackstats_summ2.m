%hlid_mds_transform_jackstats_summ2: summary plots for hlid_mds_transform_jackstats
% focusing on directions of principal axes, and max Euclidean dimension, jackknifed on stimuli
% (embeddings via pca, mds, cosine distances, and Pearson distances)
%
%runs on results structure from hlid_mds_transform_jackstats
%
% plottingg code borrowed from hlid_geom_transform_stats_summ, adapted for
% jackknife (so no selection based on mixing entropy, and no randomizations)
%
%  See also:  HLID_SETUP, HLID_MDS_TRANSFORM_JACKSTATS, HLID_GEOM_TRANSFORM_STATS_SUMM.,
%
rng_state=rng;
if (results.if_frozen~=0) 
    rng('default');
end
colors_jack=rand(results.nstims,3);
rng(rng_state);
%for compatibilty with hlid_geom_transform_stats_summ (plot_mode=1 for outputs of hlid_geom_transform_stats, 2 for hlida_mds_transform_stats)
%getinp('1 for output of hlid_geom_transform_stats (nsubs, npreprocs), 2 for hlid_mds_transform_stats (isubmean, nmeths)','d',[ 1 2]);
plot_mode=2; 
nu1=length(results.submean_use_list);
nu1_labels_all={'sm=0','sm=1'};
nu1_labels=nu1_labels_all(1+results.submean_use_list);
nu2=results.nmeths;
nu2_labels=results.meth_names_short;
if_smallfigs=getinp('1 for smaller figs','d',[0 1]); 
%
for igp=1:results.ngps
    file_range=[1:results.nsets_gp(igp)]+sum(results.nsets_gp(1:igp-1));
    med=squeeze(min(results.max_euc_dim(file_range,:,:),[],1));
    med_jack=squeeze(min(min(results.max_euc_dim_jackstim,[],1),[],4));
    disp(sprintf(' group %1.0f: max Euclidean dimension (minimum across files)',igp))
    disp('                           data     min(jackknife by stim)')
    disp('                       sm=0  sm=1       sm=0  sm=1')
    for imeth_ptr=1:length(results.meth_use_list)
        imeth=results.meth_use_list(imeth_ptr);
        disp(sprintf('%20s  %4.0f  %4.0f       %4.0f  %4.0f',results.meth_names_short{imeth},med(imeth,:),med_jack(imeth,:)));
    end
end
%
% hlid_geom_transform_stats_summ: summary plots from hlid_geom_transform_stats
%
if ~exist('if_flip_projs') if_flip_projs=1; end
%plotting params
hlid_setup;
if ~exist('display_order_spec')
    display_order=display_orders.kcmerge;
elseif isempty(display_order_spec)
    display_order=results.stimulus_names_display;
elseif ~iscell(display_order_spec)
    display_order=display_orders.(display_order_spec);
end
if ~exist('ebw') ebw=0.1; end %error bar width
if ~exist('colors') colors={'k','b','m','r',[.75 .75 0],[0 .5 0]}; end
if ~exist('hoff_withinproj') hoff_withinproj=0.12; end %horiz offset for each dimension within the rows of the projection plot
if ~exist('voff_withinproj') voff_withinproj=.14; end %vertical offset for each dimension in histogram plots
if ~exist('proj_lw') proj_lw=1; end
if ~exist('hist_bins') hist_bins=25; end
if ~exist('magnif_hist_range') magnif_hist_range=[.15 1.5]; end
%
ra_text={'ref','adj'};
%
results.if_flip_projs=if_flip_projs;
%
rbase=results.geo{1,1,1}{results.dimlist(1),results.dimlist(1)};
if (plot_mode==2)
    rbase.ref_file=rbase.ref_file(1:2+strfind(rbase.ref_file,'...'));
    rbase.adj_file=rbase.adj_file(1:2+strfind(rbase.adj_file,'...'));
end
dimlist=results.dimlist;
model_types=rbase.model_types_def.model_types;
%
%plot axis magnifications for adj dim = ref dim
results.magnif_summ_dims='d1: nsubs or ifsubmean, d2: npreprocs or nmeths, d3: nembeds, d4: nmodels';
results.magnif_summ_dims_inside='d1: ref dim= adj dim, [d2: eiv number] [boot or shuffle]';
results.magnif_summ=cell(nu1,nu2,results.nembeds,results.nmodels);
%
results.projs_summ_dims='d1: nsubs or if submean, d2: npreprocs or nments, d3: nembeds, d4: nmodels';
results.projs_summ_dims_inside='d1: ref dim= adj dim, d2: 1->ref, 2->adj; then stim and coord and [boot]';
results.projs_summ=cell(nu1,nu2,results.nembeds,results.nmodels);
%
for imodel=1:results.nmodels
    for iembed=1:results.nembeds
        for iu1=1:nu1
            for iu2=1:nu2
                %
                %collect for each model dimension (ref dim = adj dim = idim), then do vectorized calcs
                %for jackknifes, collect magnif factors and projections, keeping projections on ref and also on adj
                %
                magnif_all=NaN(max(dimlist),max(dimlist));
                magnif_all_jack=NaN(max(dimlist),max(dimlist),results.nstims);
                %
                projs=cell(max(dimlist),2); %dim 2: 1 is ref, 2 is adj: 
                projs_jack=cell(max(dimlist),2);
                for idim_ptr=1:length(dimlist)
                    idim=dimlist(idim_ptr);
                    %get magnif factors for full data and jackknifes
                    magnif_all(dimlist(idim_ptr),[1:idim])=results.geo_majaxes{iu1,iu2,iembed}{idim,idim}.ref.magnifs{imodel}';
                    %
                    for ijack=1:results.nstims
                        magnif_all_jack(dimlist(idim_ptr),[1:idim],ijack)=results.geo_majaxes_jack_by_stim{iu1,iu2,iembed,ijack}{idim,idim}.ref.magnifs{imodel}';
                    end %ijack
                    %get projections for full data and jackknifes
                    for ira=1:2
                        projs_orig=results.geo_majaxes{iu1,iu2,iembed}{idim,idim}.(ra_text{ira}).projections{imodel};
                        projs{idim,ira}=projs_orig;
                        if if_flip_projs
                            for ieiv=1:idim
                                flip=sign(projs_orig(find(abs(projs_orig(:,ieiv))==max(abs(projs_orig(:,ieiv)))),ieiv)); %largest abs value made positive
                                projs{idim,ira}(:,ieiv)=projs_orig(:,ieiv)*flip;
                            end
                        end
                        projs_jack{idim,ira}=NaN(size(projs{idim,ira},1),idim,results.nstims);
                        pjack=NaN(size(projs{idim,ira},1),idim,2);
                    end
                    for ijack=1:results.nstims
                        keep=setdiff([1:results.nstims],ijack);
                        jack_base=results.geo_majaxes_jack_by_stim{iu1,iu2,iembed,ijack}{idim,idim};
                %         %flip the bootstrapped projections if needed to align with non-bootstrapped data
                        for ira=1:2
                            pjack(keep,:,ira)=jack_base.(ra_text{ira}).projections{imodel};
                            pjack(ijack,:,ira)=NaN;
                            %these two lines are a key fix
                            flips=sign(sum(pjack(keep,:,ira).*projs{idim,ira}(keep,:),1)); %dot-product omitting the jackknifed stimuls
                            projs_jack{idim,ira}(:,:,ijack)=repmat(flips,results.nstims,1).*pjack(:,:,ira);
                        end %ira
                    end %ijack
                end %idim_ptr
                %
                m=struct;
                m.magnif_all=magnif_all;
                m.magnif_all_jack=magnif_all_jack;
                results.magnif_summ{iu1,iu2,iembed,imodel}=m;
                %
                p=struct;
                p.projs=projs;
                p.projs_jack=projs_jack;
                results.projs_summ{iu1,iu2,iembed,imodel}=p;
            end %iu2
        end %iu1
    end %iembed
end %imodel
for imodel=1:results.nmodels
    for iembed=1:results.nembeds
        %
        %projection plots
        %
        % for ira=1:2
        %     figure;
        %     if if_smallfigs
        %         set(gcf,'Position',[100 100 1200 800]);
        %     else
        %         set(gcf,'Position',[100 50 1600 950]);
        %     end
        %     set(gcf,'NumberTitle','off');
        %     set(gcf,'Name',cat(2,'major axes for ',ra_text{ira},' ',model_types{imodel},', ',results.embed_labels{iembed}));
        %     for iu1=1:nu1
        %         for iu2=1:nu2
        %             m=results.magnif_summ{iu1,iu2,iembed,imodel};
        %             p=results.projs_summ{iu1,iu2,iembed,imodel};
        %             subplot(nu1,nu2,iu1+(iu2-1)*nu1);
        %             for idim_ptr=1:length(dimlist)
        %                 idim=dimlist(idim_ptr);
        %                 yplot_off=2*idim-1;
        %                 plot([0.5 results.nstims+0.5],repmat(yplot_off,1,2),'k:');
        %                 hold on;
        %                 %plot the projections
        %                 proj_plot=p.projs{idim,ira};
        %                 for dproj=1:idim
        %                     proj_max=max(abs(proj_plot(:,dproj))); %rescale by maximum projection size
        %                     color=colors{mod(dproj-1,length(colors))+1};
        %                     xplot_off=hoff_withinproj*(dproj-0.5*(idim+1));
        %                     for istim_name=1:length(display_order) %stim_name determines position to plot
        %                         istim=strmatch(display_order{istim_name},results.stimulus_names_display,'exact');
        %                         if ~isempty(istim)
        %                             plot(repmat(xplot_off+istim_name,2,1),yplot_off+[0 proj_plot(istim,dproj)/proj_max],'Color',color,'LineWidth',proj_lw);
        %                             for ib=1:nbq
        %                                 plot(repmat(xplot_off+istim_name,1,2)+hoff_withinproj*[-0.5 0.5],...
        %                                     yplot_off+repmat(quantile(squeeze(p.projs_boot{idim,ira}(istim,dproj,:)),boot_quantiles(ib))/proj_max,1,2),...
        %                                     'Color',color,'LineWidth',1);
        %                             end %ibq
        %                         end %isempty
        %                     end %stim_name
        %                 end %dproj (projection dimension)
        %             end %idim_ptr
        %             set(gca,'XLim',[1 results.nstims]+[-0.5 0.5]);
        %             set(gca,'XTick',[1:results.nstims]);
        %             set(gca,'XTickLabel',display_order);
        %             set(gca,'YLim',[0 2*max(dimlist)]);
        %             set(gca,'YTick',[1:2:2*max(dimlist)-1]);
        %             set(gca,'YTickLabel',[1:max(dimlist)]);
        %             ylabel(sprintf('%s dim',ra_text{ira}));
        %             title(sprintf('%s %s',nu1_labels{iu1},nu2_labels{iu2}),'Interpreter','none')
        %         end %iu2
        %     end %iu1
        %     hlid_geom_transform_stats_label;
        %     axes('Position',[0.75,0.05,0.01,0.01]);
        %     text(0,0,sprintf('flip projections so max is >0: %1.0f',if_flip_projs));
        %     axis off
        % end %ira
        % %
        % %plots of magnif ratio distributions and cosine distributions
        % %
        % nsp=3; %3 kinds of subplots: dot prod in ref space, dot prod in adj space, magnif factor 
        % if (plot_mode==1)           
        %     nsp_together=3;
        % else
        %     nsp_together=1;
        % end
        % for isp_lo=1:3/nsp_together
        %     figure;
        %     if if_smallfigs
        %         set(gcf,'Position',[100 100 1200 800]);
        %     else
        %         set(gcf,'Position',[100 50 1600 950]);
        %     end
        %     set(gcf,'NumberTitle','off');
        %     set(gcf,'Name',cat(2,'ratio and cosine distribs for ',ra_text{ira},' ',model_types{imodel},', ',results.embed_labels{iembed}));
        %     for iu1=1:nu1
        %         for iu2=1:nu2
        %             m=results.magnif_summ{iu1,iu2,iembed,imodel};
        %             p=results.projs_summ{iu1,iu2,iembed,imodel};
        %             %three plots: dot prod in ref space, dot prod in adj space, magnif ratio
        %             for isp=isp_lo:isp_lo+nsp_together-1;
        %                 subplot(nu1,nu2*nsp_together,(isp-isp_lo+1)+nsp_together*(iu1+(iu2-1)*nu1-1));
        %                 for idim_ptr=1:length(dimlist)
        %                     idim=dimlist(idim_ptr);
        %                     yplot_off=idim-1;
        %                     switch isp
        %                         case {1,2} %dot products
        %                             xlims=[0 1];
        %                             if_semilogx=0;
        %                             xlabel_string=cat(2,'cosine ',ra_text{isp});
        %                             bin_edges_plot=[0:hist_bins]*(1/hist_bins);
        %                         case 3 %magnif factor
        %                             xlims=magnif_hist_range;
        %                             xlabel_string='magnif factor';
        %                             if_semilogx=1;
        %                             bin_edges_plot=magnif_hist_range(1)*(magnif_hist_range(2)/magnif_hist_range(1)).^([0:hist_bins]/hist_bins);
        %                     end
        %                     bin_edges=[-Inf bin_edges_plot(1:end-1) Inf];
        %                     plot(xlims,repmat(yplot_off,1,2),'k:');
        %                     hold on;
        %                     for dproj=1:idim
        %                         color=colors{mod(dproj-1,length(colors))+1};
        %                         switch isp
        %                             case {1,2} %norm dot products
        %                                 proj_exp=p.projs{idim,isp}(:,dproj);
        %                                 proj_boot=reshape(p.projs_boot{idim,isp}(:,dproj,:),size(p.projs_boot{idim,isp},1),results.nboots_within);
        %                                 data_exp=1;
        %                                 data_boot=(proj_exp'*proj_boot)./sqrt(sum(proj_exp.^2))./sqrt(sum(proj_boot.^2,1));
        %                             case 3 %magnif factor
        %                                 data_exp=m.magnif_all(idim,dproj); %01Jan26: first index was wrongly dproj
        %                                 data_boot=squeeze(m.magnif_all_boot(idim,dproj,:)); %01Jan26: first index was wrongly dproj
        %                         end
        %                         %plot histogram "by hand", so we can control offset
        %                         hist_data_raw=histc(data_boot,bin_edges)';
        %                         hist_data=zeros(hist_bins,1);
        %                         hist_data(1)=hist_data_raw(1)+hist_data_raw(2);
        %                         hist_data(2:hist_bins)=hist_data_raw(3:end-1);
        %                         hist_data=hist_data/max(hist_data(:));
        %                         hist_xvals=reshape(repmat(bin_edges_plot,2,1),2*hist_bins+2,1);
        %                         hist_yvals=reshape(repmat(hist_data(:)',2,1),2*hist_bins,1);
        %                         hist_yvals=[0;hist_yvals;0];
        %                         plot(hist_xvals,yplot_off+hist_yvals*0.8,'Color',color,'LineWidth',1);
        %                         plot(data_exp,yplot_off+voff_withinproj*dproj,'*','Color',color); %plot actual value from data
        %                         plot(quantile(data_boot,boot_quantiles([1 end])),repmat(yplot_off+voff_withinproj*dproj,1,2),'Color',color); %show the inter-quantile range
        %                         for ibq=1:nbq
        %                             plot(repmat(quantile(data_boot,boot_quantiles(ibq)),1,2),yplot_off+voff_withinproj*(dproj+[-0.5 0.5]),'Color',color);
        %                         end
        %                     end %dproj (projection dimension)
        %                 end %idim_ptr
        %                 set(gca,'YLim',[0 max(dimlist)]);
        %                 set(gca,'YTick',[0:max(dimlist)-1]);
        %                 set(gca,'YTickLabel',[1:max(dimlist)]);
        %                 ylabel('dim');
        %                 if if_semilogx
        %                     set(gca,'XScale','log');
        %                     set(gca,'XTick',unique([xlims 1]));
        %                 else
        %                     set(gca,'XTick',unique([xlims mean(xlims)]));
        %                 end
        %                 set(gca,'XLim',xlims);
        %                 xlabel(xlabel_string);
        %                 if plot_mode==1
        %                     if isp==2
        %                         title(sprintf('%s %s',nu1_labels{iu1},nu2_labels{iu2}),'Interpreter','none');
        %                     end
        %                 else
        %                     title(sprintf('%s %s',nu1_labels{iu1},nu2_labels{iu2}),'Interpreter','none');
        %                 end
        %             end %isp
        %         end %iu2
        %     end %iu1
        %     hlid_geom_transform_stats_label;
        %     axes('Position',[0.75,0.05,0.01,0.01]);
        %     text(0,0,sprintf('flip projections so max is >0: %1.0f',if_flip_projs));
        %     axis off
        % end %isp_lo
    end %iembed
end %imodel
%

