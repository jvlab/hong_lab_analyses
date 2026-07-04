%zheng_apl_align_stats_rs: statistics of alignment, TNT and control files
% Modular version of zheng_apl_align_stats
%run after zheng_apl_read_align_rs.m, adapted from psg_align_stats_demo.m
%
% See also: RS_KNIT_COORDSETS.
%
if_debug=getinp('1 for debug mode','d',[0 1],0);
%
if_savemat=getinp('1 to save key variables in mat-file','d',[0 1]);
if_savefig=getinp('1 to save figure(s)','d',[0 1]);
file_name_base=cat(2,'zheng_apl_align_stats_',data_use,sm_string);
if (if_savemat | if_savefig)
    if_ok=0;
    while (if_ok==0)
        file_name_base=getinp('file name base','s',[],file_name_base);
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
%
%
if ~exist('shuff_quantiles') shuff_quantiles=[0.01 0.05 0.5 0.95 0.99]; end %quantiles for showing shuffled data
nquantiles=length(shuff_quantiles);
%
if ~exist('label_shorten') label_shorten={'coords','hlid_','odor17_','megamat0','6pt','3pt','sess01_10','sess01_20','__','_'}; end %strings to remove from labels
if ~exist('label_replace') label_replace={''      ,''     ,''       ,''        ,''   ,''   ,''         ,''         ,'_' ,'-'}; end %strings to replace
%
if if_debug==1
    if ~exist('nshuffs') nshuffs=10; end
else
    if ~exist('nshuffs') nshuffs=500; end
end
pcon_dim_max=nstims_all-dim_reduce;
dim_list_in=[1:pcon_dim_max];
dim_list_out=dim_list_in;
dim_list_in_max=pcon_dim_max;
dim_list_out_max=pcon_dim_max;
%
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.max_iters=1000;
opts_knit.pcon_init_method=0;
opts_knit.nshuffs=nshuffs;
opts_knit.if_frozen=1;
opts_knit.if_log=1;
%
%make short labels
%
dataset_labels=cell(1,nsets);
for iset=1:nsets
    dataset_labels{iset}=data_aligned.sets{iset}.label(1+max([find(data_aligned.sets{iset}.label=='/'),find(data_aligned.sets{iset}.label=='\')]):end);
    if exist('label_shorten')
        for id=1:length(label_shorten)
            dataset_labels{iset}=strrep(dataset_labels{iset},label_shorten{id},label_replace{id});
        end
    end
    if dataset_labels{iset}(end)=='-'
        dataset_labels{iset}=dataset_labels{iset}(1:end-1);
    end
end
%
ra_setup=struct; %use as starting point for results and for plotting
%
ra_setup.dataset_labels=dataset_labels;
ra_setup.stimulus_labels=data_aligned.sas{1}.typenames;
ra_setup.sets=data_aligned.sets;
ra_setup.data_orig=data_aligned.ds;
ra_setup.metadata_orig=data_read.sas;
ra_setup.sa_consensus=data_aligned.sas{1}; %metadata for ds_consensus
ra_setup.sas_components=data_aligned.sas;
%
ra_setup.nstims=nstims_all;
ra_setup.nsets=nsets;
ra_setup.nshuffs=nshuffs;
ra_setup.shuff_quantiles=shuff_quantiles; 
%
ra_setup.dim_list_in_max=dim_list_in_max;
ra_setup.dim_list_out_max=dim_list_out_max;
%
ra=cell(1,2);
fig_handle=figure;
set(gcf,'Position',[100 100 1300 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',cat(2,'consensus analysis: ',data_use,sm_string));
%
for allow_scale=0:1
    ia=allow_scale+1;
    %
    opts_knit.if_stats=1;
    opts_knit.if_plot=0;
    opts_knit.allow_scale=allow_scale;
    %
    aux=struct;
    aux.opts_knit=opts_knit;
    %
    [data_knit_out,aux_knit_out]=rs_knit_coordsets(data_aligned,aux);
    %save and replot
    ra{ia}=aux_knit_out.knit_stats;
    aux_replot=aux;
    aux_replot.knit_stats=aux_knit_out.knit_stats;
    aux_replot.knit_stats_setup=aux_knit_out.knit_stats_setup;
    aux_replot.knit_stats_setup.fig_handle=fig_handle;
    aux_replot.knit_stats_setup.nrows=2;
    aux_replot.knit_stats_setup.row=ia;
    aux_replot.knit_stats_setup.dataset_labels=dataset_labels;
    rs_knit_coordsets(data_aligned,aux_replot);
end
%plots; display rms dev unexplained, p-values, and quantiles
for allow_scale=0:1
    ia=allow_scale+1;
    %
    disp(sprintf('rms deviation analysis, allow_scale=%1.0f, nshuffs=%5.0f',allow_scale,ra_setup.nshuffs));
    disp('     dim  rms_avail rms unexp(data)  p   ms unexp (shuffle last dimension, quantiles)');
    quantile_text=sprintf('%8.4f  ',ra_setup.shuff_quantiles);
    disp(sprintf('%s',cat(2,repmat(' ',1,42),quantile_text)));
    %dim5 of rmsdev_overall_shuff==1 to select shuffle only last dimension
    quantiles=reshape(quantile(ra{ia}.rmsdev_overall_shuff(:,1,1,:,1),ra_setup.shuff_quantiles,4),[length(dim_list_in),length(ra_setup.shuff_quantiles)]);
    for k=1:length(dim_list_in)
        p_val=sum(double(ra{ia}.rmsdev_overall(k)>ra{ia}.rmsdev_overall_shuff(k,1,1,:,1))/ra_setup.nshuffs);
        vals=cat(2,k,ra{ia}.rmsavail_overall(k,1),ra{ia}.rmsdev_overall(k),p_val,quantiles(k,:));
        disp(vals)
    end
end
%save files
if (if_savemat)
    file_name_mat=file_name_base;
    save(file_name_mat,'ra','ra_setup','opts*','-v7.3');
    disp(sprintf('variables saved as %s',file_name_mat));
end
if (if_savefig)
    file_name_fig=file_name_base;
    savefig(gcf,file_name_fig);
    disp(sprintf('figure saved as %s',file_name_fig));
end
