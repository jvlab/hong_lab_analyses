%hlid_geom_transform_stats_label: utility labeling for hlid_geom_transform_stats
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
text(0,0,cat(2,'shuffle quantiles:',sprintf(' %6.4f',ratio_quantiles),sprintf('; bootstrap error bars (%4.0f boots):',results.nboots_within),sprintf(' %6.4f',boot_quantiles)));
axis off
axes('Position',[0.5,0.01,0.01,0.01]);
text(0,0,sprintf('shuffles selected: %5.0f of %5.0f (frac: %5.3f), nrelabeled: %2.0f to %2.0f',...
    length(mixents_ptrs),length(mixents),mixent_frac,min(nrelabeled(mixents_ptrs)),max(nrelabeled(mixents_ptrs))));
axis off;
