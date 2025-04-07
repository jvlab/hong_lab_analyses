%hlid_geom_transform_stats_summ: summary plots from hlid_geom_transform_statslength(dimlist)
% 
% runs on results structure from hlid_geom_transform_stats
% to plot and tabulate major axis ratios, for adj dim = ref dim
%
% projections are plotted after normalization so that max projection = 1
% and are optoinally flipped (by if_flip_projs) so that largest projection is positive
%
% ratio_quantiles can be customized for significance levels
% if_show_max can be set to 1 to also show max axis ratio
% mixent_frac can be set to less than 1 to only inclucde shuffles with a minimum amount of mixing entropies
% if_flip_projs can be set to 0 to disable flipping of projections so that largest projection is positive
%   (boostrapped projections are aligned b dot-product to unbootstrapped versions)
% display_order_spec can be used to specify the order of stimuli in the display
%   defaults to display_orders.kcmerge from hlid_setup
%   display_order_spec=[] to for native order
%   display_oreer_spec='kcmerge','kctnt', or 'kcclust' to use one of the lists in display_orders, from hlid_setup.
%
%   See also:  HLID_GEOM_TRANSFORM_STATS, MULII_SHUFF_MIXENT,
%   HLID_GEOM_TRANSFORM_STATS_LABEL.
%
if ~exist('ratio_quantiles') ratio_quantiles=[.5 .05 .01]; end %to look at top of distribution
nrq=length(ratio_quantiles);
if ~exist('if_show_max') if_show_max=0; end % set to 1 to show maximum axis ratio
if ~exist('mixent_frac') mixent_frac=1; end %set to < 1 to only include shuffles with a large mixing entropy
if ~exist('boot_quantiles') boot_quantiles=[0.025 0.975]; end %make empty to omit bootstraps
nbq=length(boot_quantiles);
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
rng_state=rng;
if (results.if_frozen~=0) 
    rng('default');
end
%
[mixents,mixmats]=multi_shuff_mixent(results.shuffs_between,results.gps);
nrelabeled=zeros(results.nshuffs_between,1);
for igp=1:results.ngps
    nrelabeled=nrelabeled+sum(results.shuff_gp_origs{igp}~=igp,2);
end
mixents_maxposs=-sum(lognz(results.nsets_gp/results.nsets).*(results.nsets_gp/results.nsets))/log(2);
disp(sprintf('before selection: %5.0f shuffles, mixing entropies range from %7.4f to %7.4f, %4.0f unique values, max possible is %7.4f, nrelabeled: %2.0f to %2.0f',...
    length(mixents),min(mixents),max(mixents),length(unique(mixents)),mixents_maxposs,min(nrelabeled),max(nrelabeled)));
%add a small amount of random jitter so that quantiles will give close to the exact number desired
mixent_jit=min([1;abs(diff(unique(mixents)))]);
mixents_jittered=mixents+rand(results.nshuffs_between,1)*mixent_jit;
mixents_minkeep=quantile(mixents_jittered,[1-mixent_frac]);
mixents_ptrs=find(mixents_jittered>=mixents_minkeep);
%
disp(sprintf(' after selection: %5.0f shuffles, mixing entropies range from %7.4f to %7.4f, %4.0f unique values, max possible is %7.4f, nrelabeled: %2.0f to %2.0f',...
    length(mixents(mixents_ptrs)),min(mixents(mixents_ptrs)),max(mixents(mixents_ptrs)),length(unique(mixents(mixents_ptrs))),mixents_maxposs, ...
    min(nrelabeled(mixents_ptrs)),max(nrelabeled(mixents_ptrs))));
%
if ~isfield(results,'nboots_within') %for compatibility
    results.nboots_within=0;
end
%
rbase=results.geo{1,1,1}{results.dimlist(1),results.dimlist(1)};
dimlist=results.dimlist;
model_types=rbase.model_types_def.model_types;
%
%plot axis magnifications for adj dim = ref dim
results.magnif_summ_dims='d1: nsubs, d2: npreprocs, d3: nembeds, d4: nmodels';
results.magnif_summ_dims_inside='d1: ref dim= adj dim, [d2: eiv number] [boot or shuffle]';
results.magnif_summ=cell(results.nsubs,results.npreprocs,results.nembeds,results.nmodels);
%
results.projs_summ_dims='d1: nsubs, d2: npreprocs, d3: nembeds, d4: nmodels';
results.projs_summ_dims_inside='d1: ref dim= adj dim, d2: 1->ref, 2->adj; then stim and coord and [boot]';
results.projs_summ=cell(results.nsubs,results.npreprocs,results.nembeds,results.nmodels);
%
for imodel=1:results.nmodels
    for iembed=1:results.nembeds
        for isub=1:results.nsubs
            for ipreproc=1:results.npreprocs
                %
                %collect for each model dimension (ref dim = adj dim = idim), then do vectorized calcs
                %for shuffles, only collect magnif factors
                %for boots, collect magnif factors and projections, keeping projections on ref and also on adj
                %
                magnif_all=NaN(max(dimlist),max(dimlist));
                magnif_rng=NaN(max(dimlist),1);
                %
                magnif_all_shuff=NaN(max(dimlist),max(dimlist),results.nshuffs_between);
                magnif_rng_shuff=NaN(max(dimlist),results.nshuffs_between);
                %
                magnif_all_boot=NaN(max(dimlist),max(dimlist),results.nboots_within);
                magnif_rng_boot=NaN(max(dimlist),results.nboots_within);
                %
                projs=cell(max(dimlist),2); %dim 2: 1 is ref, 2 is adj: 
                projs_boot=cell(max(dimlist),2);
                for idim_ptr=1:length(dimlist)
                    idim=dimlist(idim_ptr);
                    magnif_all(dimlist(idim_ptr),[1:idim])=results.geo_majaxes{isub,ipreproc,iembed}{idim,idim}.ref.magnifs{imodel}';
                    magnif_rng(dimlist(idim_ptr),1)=results.geo_majaxes{isub,ipreproc,iembed}{idim,idim}.ref.magnif_ratio{imodel};
                    %
                    for ishuff=1:results.nshuffs_between
                        magnif_all_shuff(dimlist(idim_ptr),[1:idim],ishuff)=results.geo_majaxes_shuff{isub,ipreproc,iembed,ishuff}{idim,idim}.ref.magnifs{imodel}';
                        magnif_rng_shuff(dimlist(idim_ptr),ishuff)=results.geo_majaxes_shuff{isub,ipreproc,iembed,ishuff}{idim,idim}.ref.magnif_ratio{imodel};
                    end %ishuff
                    %
                    for ira=1:2
                        projs_orig=results.geo_majaxes{isub,ipreproc,iembed}{idim,idim}.(ra_text{ira}).projections{imodel};
                        projs{idim,ira}=projs_orig;
                        if if_flip_projs
                            for ieiv=1:idim
                                flip=sign(projs_orig(find(abs(projs_orig(:,ieiv))==max(abs(projs_orig(:,ieiv)))),ieiv));
                                projs{idim,ira}(:,ieiv)=projs_orig(:,ieiv)*flip;
                            end
                        end
                        projs_boot{idim,ira}=NaN(size(projs{idim,ira},1),idim,results.nboots_within);
                        pboot=NaN(size(projs{idim,ira},1),idim,2);
                    end
                    for iboot=1:results.nboots_within
                        boot_base=results.geo_majaxes_boot{isub,ipreproc,iembed,iboot}{idim,idim};
                        magnif_all_boot(dimlist(idim_ptr),[1:idim],iboot)=boot_base.ref.magnifs{imodel}';
                        magnif_rng_boot(dimlist(idim_ptr),iboot)=boot_base.ref.magnif_ratio{imodel};
                        %flip the bootstrapped projections if needed to align with non-bootstrapped data
                        for ira=1:2
                            pboot(:,:,ira)=boot_base.(ra_text{ira}).projections{imodel};
                            flips=sign(pboot(:,:,ira).*projs{idim,ira});
                            projs_boot{idim,ira}(:,:,iboot)=flips.*pboot(:,:,ira);
                        end %ira
                    end %iboot
                end %idim_ptr
                dm=min(2,size(magnif_all,2)); %in case only one dim
                magnif_r12=magnif_all(:,1)./magnif_all(:,dm);
                magnif_rgm=magnif_all(:,1)./geomean(magnif_all,2,'omitnan');
                magnif_r12_shuff=reshape(magnif_all_shuff(:,1,:)./magnif_all_shuff(:,dm,:),[max(dimlist) results.nshuffs_between]);
                magnif_rgm_shuff=reshape(magnif_all_shuff(:,1,:)./geomean(magnif_all_shuff,2,'omitnan'),[max(dimlist) results.nshuffs_between]);
                magnif_r12_boot=reshape(magnif_all_boot(:,1,:)./magnif_all_boot(:,dm,:),[max(dimlist) results.nboots_within]);
                magnif_rgm_boot=reshape(magnif_all_boot(:,1,:)./geomean(magnif_all_boot,2,'omitnan'),[max(dimlist) results.nboots_within]);
                %
                m=struct;
                m.magnif_all=magnif_all;
                m.magnif_rng=magnif_rng;
                m.magnif_r12=magnif_r12;
                m.magnif_rgm=magnif_rgm;
                m.magnif_all_shuff=magnif_all_shuff;
                m.magnif_rng_shuff=magnif_rng_shuff;
                m.magnif_r12_shuff=magnif_r12_shuff;
                m.magnif_rgm_shuff=magnif_rgm_shuff;
                m.magnif_all_boot=magnif_all_boot;
                m.magnif_rng_boot=magnif_rng_boot;
                m.magnif_r12_boot=magnif_r12_boot;
                m.magnif_rgm_boot=magnif_rgm_boot;
                results.magnif_summ{isub,ipreproc,iembed,imodel}=m;
                %
                p=struct;
                p.projs=projs;
                p.projs_boot=projs_boot;
                results.projs_summ{isub,ipreproc,iembed,imodel}=p;
            end %ipreproc
        end %isub
    end %iembed
end %imodel
for imodel=1:results.nmodels
    for iembed=1:results.nembeds
        %
        %magnification factor plots
        %
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'magnif factors for ',model_types{imodel},', ',results.embed_labels{iembed}));
        for isub=1:results.nsubs
            for ipreproc=1:results.npreprocs
                m=results.magnif_summ{isub,ipreproc,iembed,imodel};
                subplot(results.nsubs,results.npreprocs,isub+(ipreproc-1)*results.nsubs);               
                hl=cell(0);
                ht=[];
                hp=semilogy(dimlist,m.magnif_rng(dimlist),'k','LineWidth',2);
                hold on;
                hl=[hl,hp];
                ht=strvcat(ht,'high/low');
                if nrq>0
                    for iq=1:nrq
                        hp=semilogy(dimlist,quantile(m.magnif_rng_shuff(dimlist,mixents_ptrs),1-ratio_quantiles(iq),2),'k','LineWidth',1);
                    end
                end
                %
                if length(dimlist)>=2
                    hp=semilogy(dimlist,m.magnif_r12(dimlist),'r','LineWidth',2);
                    hl=[hl,hp];
                    ht=strvcat(ht,'high/second');
                    hp=semilogy(dimlist,m.magnif_rgm(dimlist),'g','LineWidth',2);
                    hl=[hl,hp];
                    ht=strvcat(ht,'high/geomean');
                    if if_show_max  
                        hp=semilogy(dimlist,m.magnif_all(dimlist,1),'m','LineWidth',2);
                        hl=[hl,hp];
                        ht=strvcat(ht,'high');
                    end
                    if nrq>0
                        for iq=1:nrq
                            hp=semilogy(dimlist,quantile(m.magnif_r12_shuff(dimlist,mixents_ptrs),1-ratio_quantiles(iq),2),'r','LineWidth',1);
                            hp=semilogy(dimlist,quantile(m.magnif_rgm_shuff(dimlist,mixents_ptrs),1-ratio_quantiles(iq),2),'g','LineWidth',1);
                            if if_show_max                                                               
                                hp=semilogy(dimlist,quantile(m.magnif_all_shuff(dimlist,1,mixents_ptrs),1-ratio_quantiles(iq),3),'m','LineWidth',1);
                                hp=semilogy(dimlist,quantile(m.magnif_all_shuff(dimlist,1,mixents_ptrs),ratio_quantiles(iq),3),'m:','LineWidth',1);
                            end
                        end %iq
                    end %nrq
                    if nbq>1 & results.nboots_within>0
                        boot_ebs_rng=quantile(m.magnif_rng_boot,boot_quantiles,2);
                        for idim_ptr=1:length(dimlist)
                            idim=dimlist(idim_ptr);
                            semilogy(repmat(idim,1,nbq),boot_ebs_rng(idim,:),'k');
                            for ib=1:nbq
                                semilogy([idim idim]+ebw*[-1 1],repmat(boot_ebs_rng(idim,ib),1,2),'k')
                            end
                        end
                        if length(dimlist)>=2
                            boot_ebs_r12=quantile(m.magnif_r12_boot,boot_quantiles,2);
                            boot_ebs_rgm=quantile(m.magnif_rgm_boot,boot_quantiles,2);
                            boot_ebs_all=quantile(squeeze(m.magnif_all_boot(:,1,:)),boot_quantiles,2);
                            for idim_ptr=1:length(dimlist)
                                idim=dimlist(idim_ptr);
                                semilogy(repmat(idim,1,nbq),boot_ebs_r12(idim,:),'r');
                                semilogy(repmat(idim,1,nbq),boot_ebs_rgm(idim,:),'g');
                                for ib=1:nbq
                                    semilogy([idim idim]+ebw*[-1 1],repmat(boot_ebs_r12(idim,ib),1,2),'r');
                                    semilogy([idim idim]+ebw*[-1 1],repmat(boot_ebs_rgm(idim,ib),1,2),'g')
                                end
                                if if_show_max
                                    semilogy(repmat(idim,1,nbq),boot_ebs_all(idim,:),'m');
                                    for ib=1:nbq
                                        semilogy([idim idim]+ebw*[-1 1],repmat(boot_ebs_all(idim,ib),1,2),'m');
                                    end
                                end %if_show_max
                            end %idim_ptr
                        end %dimlist
                    end %have bootstraps?
                end
                %
                legend(hl,ht,'Location','NorthWest');
                set(gca,'XLim',[-0.5 0.5]+[1 max(dimlist)]);
                set(gca,'XTick',[1:max(dimlist)]);
                set(gca,'YLim',[1 10]);
                xlabel('dimension');
                ylabel('ratio');
                title(sprintf('%s %s',results.sub_labels{isub},results.preproc_labels{ipreproc}),'Interpreter','none')
             end %ipreproc
        end %isub
        hlid_geom_transform_stats_label;
        %
        %projection plots
        %
        for ira=1:2
            figure;
            set(gcf,'Position',[100 50 1600 950]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'major axes for ',ra_text{ira},' ',model_types{imodel},', ',results.embed_labels{iembed}));
            for isub=1:results.nsubs
                for ipreproc=1:results.npreprocs
                    m=results.magnif_summ{isub,ipreproc,iembed,imodel};
                    p=results.projs_summ{isub,ipreproc,iembed,imodel};
                    subplot(results.nsubs,results.npreprocs,isub+(ipreproc-1)*results.nsubs);
                    for idim_ptr=1:length(dimlist)
                        idim=dimlist(idim_ptr);
                        yplot_off=2*idim-1;
                        plot([0.5 results.nstims+0.5],repmat(yplot_off,1,2),'k:');
                        hold on;
                        %plot the projections
                        proj_plot=p.projs{idim,ira};
                        for dproj=1:idim
                            proj_max=max(abs(proj_plot(:,dproj))); %rescale by maximum projection size
                            color=colors{mod(dproj-1,length(colors))+1};
                            xplot_off=hoff_withinproj*(dproj-0.5*(idim+1));
                            for istim_name=1:length(display_order) %stim_name determines position to plot
                                istim=strmatch(display_order{istim_name},results.stimulus_names_display,'exact');
                                if ~isempty(istim)
                                    plot(repmat(xplot_off+istim_name,2,1),yplot_off+[0 proj_plot(istim,dproj)/proj_max],'Color',color,'LineWidth',proj_lw);
                                    for ib=1:nbq
                                        plot(repmat(xplot_off+istim_name,1,2)+hoff_withinproj*[-0.5 0.5],...
                                            yplot_off+repmat(quantile(squeeze(p.projs_boot{idim,ira}(istim,dproj,:)),boot_quantiles(ib))/proj_max,1,2),...
                                            'Color',color,'LineWidth',1);
                                    end %ibq
                                end %isempty
                            end %stim_name
                        end %dproj (projection dimension)
                    end %idim_ptr
                    set(gca,'XLim',[1 results.nstims]+[-0.5 0.5]);
                    set(gca,'XTick',[1:results.nstims]);
                    set(gca,'XTickLabel',display_order);
                    set(gca,'YLim',[0 2*max(dimlist)]);
                    set(gca,'YTick',[1:2:2*max(dimlist)-1]);
                    set(gca,'YTickLabel',[1:max(dimlist)]);
                    ylabel(sprintf('%s dim',ra_text{ira}));
                    title(sprintf('%s %s',results.sub_labels{isub},results.preproc_labels{ipreproc}),'Interpreter','none')
                end %ipreproc
            end %isub
            hlid_geom_transform_stats_label;
            axes('Position',[0.75,0.05,0.01,0.01]);
            text(0,0,sprintf('flip projections so max is >0: %1.0f',if_flip_projs));
            axis off
        end %ira
        %
        %plots of magnif ratio distributions and cosine distributions
        %
        figure;
        set(gcf,'Position',[100 50 1600 950]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'ratio and cosine distribs for ',ra_text{ira},' ',model_types{imodel},', ',results.embed_labels{iembed}));
        for isub=1:results.nsubs
            for ipreproc=1:results.npreprocs
                m=results.magnif_summ{isub,ipreproc,iembed,imodel};
                p=results.projs_summ{isub,ipreproc,iembed,imodel};
                %three plots: dot prod in ref space, dot prod in adj space, magnif ratio
                nsp=3;
                for isp=1:nsp
                    subplot(results.nsubs,results.npreprocs*nsp,isp+nsp*(isub+(ipreproc-1)*results.nsubs-1));
                    for idim_ptr=1:length(dimlist)
                        idim=dimlist(idim_ptr);
                        yplot_off=idim-1;
                        switch isp
                            case {1,2} %dot products
                                xlims=[0 1];
                                if_semilogx=0;
                                xlabel_string=cat(2,'cosine ',ra_text{isp});
                                bin_edges_plot=[0:hist_bins]*(1/hist_bins);
                            case 3 %magnif factor
                                xlims=magnif_hist_range;
                                xlabel_string='magnif factor';
                                if_semilogx=1;
                                bin_edges_plot=magnif_hist_range(1)*(magnif_hist_range(2)/magnif_hist_range(1)).^([0:hist_bins]/hist_bins);
                        end
                        bin_edges=[-Inf bin_edges_plot(1:end-1) Inf];
                        plot(xlims,repmat(yplot_off,1,2),'k:');
                        hold on;
                        for dproj=1:idim
                            color=colors{mod(dproj-1,length(colors))+1};
                            switch isp
                                case {1,2} %norm dot products
                                    proj_exp=p.projs{idim,isp}(:,dproj);
                                    proj_boot=reshape(p.projs_boot{idim,isp}(:,dproj,:),size(p.projs_boot{idim,isp},1),results.nboots_within);
                                    data_exp=1;
                                    data_boot=(proj_exp'*proj_boot)./sqrt(sum(proj_exp.^2))./sqrt(sum(proj_boot.^2,1));
                                case 3 %magnif factor
                                    data_exp=m.magnif_all(dproj,dproj);
                                    data_boot=squeeze(m.magnif_all_boot(dproj,dproj,:));
                            end
                            %plot histogram "by hand", so we can control offset
                            hist_data_raw=histc(data_boot,bin_edges)';
                            hist_data=zeros(hist_bins,1);
                            hist_data(1)=hist_data_raw(1)+hist_data_raw(2);
                            hist_data(2:hist_bins)=hist_data_raw(3:end-1);
                            hist_data=hist_data/max(hist_data(:));
                            hist_xvals=reshape(repmat(bin_edges_plot,2,1),2*hist_bins+2,1);
                            hist_yvals=reshape(repmat(hist_data(:)',2,1),2*hist_bins,1);
                            hist_yvals=[0;hist_yvals;0];
                            plot(hist_xvals,yplot_off+hist_yvals*0.8,'Color',color,'LineWidth',1);
                            plot(data_exp,yplot_off+voff_withinproj*dproj,'*','Color',color); %plot actual value from data
                            plot(quantile(data_boot,boot_quantiles([1 end])),repmat(yplot_off+voff_withinproj*dproj,1,2),'Color',color); %show the inter-quantile range
                            for ibq=1:nbq
                                plot(repmat(quantile(data_boot,boot_quantiles(ibq)),1,2),yplot_off+voff_withinproj*(dproj+[-0.5 0.5]),'Color',color);
                            end
                        end %dproj (projection dimension)
                    end %idim_ptr
                    set(gca,'YLim',[0 max(dimlist)]);
                    set(gca,'YTick',[0:max(dimlist)-1]);
                    set(gca,'YTickLabel',[1:max(dimlist)]);
                    ylabel('dim');
                    if if_semilogx
                        set(gca,'XScale','log');
                        set(gca,'XTick',unique([xlims 1]));
                    else
                        set(gca,'XTick',unique([xlims mean(xlims)]));
                    end
                    set(gca,'XLim',xlims);
                    xlabel(xlabel_string);
                    if isp==2
                        title(sprintf('%s %s %s',results.sub_labels{isub},results.preproc_labels{ipreproc}),'Interpreter','none')
                    end
                end %isp
            end %ipreproc
        end %isub
        hlid_geom_transform_stats_label;
        axes('Position',[0.75,0.05,0.01,0.01]);
        text(0,0,sprintf('flip projections so max is >0: %1.0f',if_flip_projs));
        axis off
    end %iembed
end %imodel
%
rng(rng_state); %restore random state
