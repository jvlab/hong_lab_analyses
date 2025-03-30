%hlid_geom_transform_stats_summ: summary plots from hlid_geom_transform_stats
% 
% runs on results structure from hlid_geom_transform_stats
% to plot and tabulate major axis ratios, for adj dim = ref dim
%
% ratio_quantiles can be customized for significance levelsalues
% if_show_max can be set to 1 to also show max axis ratio
% mixent_frac can be set to less than 1 to only inclucde shuffles with 
%   a minimum amount of mixing entropies
%
%   See also:  HLID_GEOM_TRANSFORM_STATS, MULII_SHUFF_MIXENT
%
if ~exist('ratio_quantiles') ratio_quantiles=[.5 .05 .01]; end %to look at top of distribution
nrq=length(ratio_quantiles);
if ~exist('if_show_max') if_show_max=0; end % set to 1 to show maximum axis ratio
if ~exist('mixent_frac') mixent_frac=1; end %set to < 1 to only include shuffles with a large mixing entropy
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
mixent_jit=min(min(abs(diff(unique(mixents)))),1);
mixents_jittered=mixents+rand(results.nshuffs_between,1)*mixent_jit;
mixents_minkeep=quantile(mixents_jittered,[1-mixent_frac]);
mixents_ptrs=find(mixents_jittered>=mixents_minkeep);
%
disp(sprintf(' after selection: %5.0f shuffles, mixing entropies range from %7.4f to %7.4f, %4.0f unique values, max possible is %7.4f, nrelabeled: %2.0f to %2.0f',...
    length(mixents(mixents_ptrs)),min(mixents(mixents_ptrs)),max(mixents(mixents_ptrs)),length(unique(mixents(mixents_ptrs))),mixents_maxposs, ...
    min(nrelabeled(mixents_ptrs)),max(nrelabeled(mixents_ptrs))));
%
rng_state=rng;
if (results.if_frozen~=0) 
    rng('default');
end
rbase=results.geo{1,1,1}{results.dimlist(1),results.dimlist(1)};
dimlist=results.dimlist;
model_types=rbase.model_types_def.model_types;
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
                        hp=semilogy(dimlist,quantile(magnif_rng_shuff(dimlist,mixents_ptrs),1-ratio_quantiles(iq),2),'k','LineWidth',1);
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
                    if if_show_max  
                        hp=semilogy(dimlist,magnif_all(dimlist,1),'m','LineWidth',2);
                        hl=[hl,hp];
                        ht=strvcat(ht,'high');
                    end
                    if nrq>0
                        magnif_r12_shuff=reshape(magnif_all_shuff(:,1,:)./magnif_all_shuff(:,2,:),[max(dimlist) results.nshuffs_between]);
                        magnif_rgm_shuff=reshape(magnif_all_shuff(:,1,:)./geomean(magnif_all_shuff,2,'omitnan'),[max(dimlist) results.nshuffs_between]);
                        for iq=1:nrq
                            hp=semilogy(dimlist,quantile(magnif_r12_shuff(dimlist,mixents_ptrs),1-ratio_quantiles(iq),2),'r','LineWidth',1);
                            hp=semilogy(dimlist,quantile(magnif_rgm_shuff(dimlist,mixents_ptrs),1-ratio_quantiles(iq),2),'g','LineWidth',1);
                            if if_show_max                                                               
                                hp=semilogy(dimlist,quantile(magnif_all_shuff(dimlist,1,mixents_ptrs),1-ratio_quantiles(iq),3),'m','LineWidth',1);
                                hp=semilogy(dimlist,quantile(magnif_all_shuff(dimlist,1,mixents_ptrs),ratio_quantiles(iq),3),'m:','LineWidth',1);
                            end
                        end %iq
                    end %nrq
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
        axes('Position',[0.01,0.03,0.01,0.01]); %for text
        text(0,0,cat(2,'ref: ',rbase.ref_file),'Interpreter','none');
        axis off       
        %
        axes('Position',[0.01,0.01,0.01,0.01]); %for text
        text(0,0,cat(2,'adj: ',rbase.adj_file),'Interpreter','none');
        axis off
        axes('Position',[0.5,0.05,0.01,0.01]);
        text(0,0,sprintf('allow scale in consensus: %1.0f',results.opts_pcon.allow_scale));
        axis off
        axes('Position',[0.5,0.03,0.01,0.01]);
        text(0,0,cat(2,'quantiles:',sprintf(' %6.4f',ratio_quantiles)));
        axis off
        axes('Position',[0.5,0.01,0.01,0.01]);
        text(0,0,sprintf('shuffles selected: %5.0f of %5.0f (frac: %5.3f), nrelabeled: %2.0f to %2.0f',...
            length(mixents_ptrs),length(mixents),mixent_frac,min(nrelabeled(mixents_ptrs)),max(nrelabeled(mixents_ptrs))));
        axis off;
    end %iembed
end %imodel
%
rng(rng_state); %restore random state
