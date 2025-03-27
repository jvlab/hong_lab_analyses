%hlid_geom_transform_stats_summ: summary plots from hlid_geom_transform_stats
% 
% runs on results structure from hlid_geom_transform_stats
% to plot and tabulate major axis ratios, for adj dim = ref dim
%
rbase=results.geo{1,1,1}{results.dimlist(1),results.dimlist(1)};
dimlist=results.dimlist;
model_types=rbase.model_types_def.model_types;
if ~exist('ratio_quantiles') ratio_quantiles=[.05 .01]; end %to look at top of distribution
nrq=length(ratio_quantiles);
%
%plot axis magnifications for adj dim = ref dim
results.magnif_summ_dims='d1: nsubs, d2: npreprocs, d3: nembeds, d4: nmodels';
results.magnif_summ_dims_inside='d1: ref dim= adj dim, d2: eiv number';
magnif_summ=cell(results.nsubs,results.npreprocs,results.nembeds,results.nmodels);
for imodel=1:results.nmodels
    for iembed=1:results.nembeds
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,model_types{imodel},', ',results.embed_labels{iembed}));
        for isub=1:results.nsubs
            for ipreproc=1:results.npreprocs
                subplot(results.nsubs,results.npreprocs,isub+(ipreproc-1)*results.nsubs);
                magnif_all=NaN(max(dimlist),max(dimlist));
                magnif_rng=NaN(max(dimlist),1);
                magnif_r12=NaN(max(dimlist),1);
                magnif_rgm=NaN(max(dimlist),1);
                %
                magnif_all_shuff=NaN(max(dimlist),max(dimlist),results.nshuffs_between);
                magnif_rng_shuff=NaN(max(dimlist),results.nshuffs_between);
                magnif_r12_shuff=NaN(max(dimlist),results.nshuffs_between);
                magnif_rgm_shuff=NaN(max(dimlist),results.nshuffs_between);
                for idim_ptr=1:length(dimlist)
                    idim=dimlist(idim_ptr);
                    magnif_all(dimlist(idim_ptr),[1:idim])=results.geo_majaxes{isub,ipreproc,iembed}{idim,idim}.ref.magnifs{imodel}';
                    magnif_rng(dimlist(idim_ptr),1)=results.geo_majaxes{isub,ipreproc,iembed}{idim,idim}.ref.magnif_ratio{imodel};
                    for ishuff=1:results.nshuffs_between
                        magnif_all_shuff(dimlist(idim_ptr),[1:idim],ishuff)=results.geo_majaxes_shuff{isub,ipreproc,iembed,ishuff}{idim,idim}.ref.magnifs{imodel}';
                        magnif_rng_shuff(dimlist(idim_ptr),ishuff)=results.geo_majaxes_shuff{isub,ipreproc,iembed,ishuff}{idim,idim}.ref.magnif_ratio{imodel};
                    end %ishuff
                end %idim_ptr
                hl=cell(0);
                ht=[];
                hp=semilogy(dimlist,magnif_rng(dimlist),'k','LineWidth',2);
                hold on;
                hl=[hl,hp];
                ht=strvcat(ht,'high/low');
                if nrq>0
                    for iq=1:nrq

                    end
                end
                %
                if length(dimlist)>=2
                    magnif_r12(:,1)=magnif_all(:,1)./magnif_all(:,2);
                    magnif_rgm(:,1)=magnif_all(:,1)./geomean(magnif_all,2,'omitnan');
                    hp=semilogy(dimlist,magnif_r12(dimlist),'r','LineWidth',2);
                    hl=[hl,hp];
                    ht=strvcat(ht,'high/second');
                    hp=semilogy(dimlist,magnif_rgm(dimlist),'g','LineWidth',2);
                    hl=[hl,hp];
                    ht=strvcat(ht,'high/geomean');
                end
                %
                legend(hl,ht,'Location','NorthWest');
                set(gca,'XLim',[-0.5 0.5]+[1 max(dimlist)]);
                set(gca,'XTick',[1:max(dimlist)]);
                set(gca,'YLim',[1 10]);
                xlabel('dimension');
                ylabel('ratio');
                title(sprintf('%s %s',results.sub_labels{isub},results.preproc_labels{ipreproc}),'Interpreter','none')
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
                results.magnif_summ{isub,ipreproc,iembed,imodel}=m;
             end %ipreproc
        end %isub
        axes('Position',[0.01,0.05,0.01,0.01]); %for text
        text(0,0,sprintf(' %s, %s, shuffles between sets: %4.0f',...
            model_types{imodel},results.embed_labels{iembed},results.nshuffs_between),'Interpreter','none');
        axis off;
        %
        axes('Position',[0.01,0.03,0.01,0.01]); %for text
        text(0,0,cat(2,'ref: ',rbase.ref_file),'Interpreter','none');
        axis off
        axes('Position',[0.01,0.01,0.01,0.01]); %for text
        text(0,0,cat(2,'adj: ',rbase.adj_file),'Interpreter','none');
        axis off
        axes('Position',[0.5,0.01,0.01,0.01]);
        text(0,0,cat(2,sprintf('allow scale: %1.0f',results.opts_pcon.allow_scale),' quantiles:',sprintf(' %6.4f',ratio_quantiles)));
        axis off
    end %iembed
end %imodel
