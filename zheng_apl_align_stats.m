%zheng_apl_align_stats: statistics of alignment, TNT and control files 
%run after zheng_apl_read_align.m, adapted from psg_align_stats_demo.m
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
if ~exist('opts_pcon') opts_pcon=struct(); end % for procrustes_consensus
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0);
opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
%
if ~exist('shuff_quantiles') shuff_quantiles=[0.01 0.05 0.5 0.95 0.99]; end %quantiles for showing shuffled data
nquantiles=length(shuff_quantiles);
%
if ~exist('label_shorten') label_shorten={'coords','hlid_','odor17_','megamat0','6pt','3pt','sess01_10','sess01_20','__','_'}; end %strings to remove from labels
if ~exist('label_replace') label_replace={''      ,''     ,''       ,''        ,''   ,''   ,''         ,''         ,'_' ,'-'}; end %strings to replace
%
if_normscale=1;
opts_pcon.if_normscale=if_normscale;
%
if if_debug==1
    if ~exist('nshuffs') nshuffs=10; end
else
    if ~exist('nshuffs') nshuffs=500; end
end
pcon_dim_max=nstims_all-dim_reduce;
dim_list_in=[1:pcon_dim_max];
dim_list_out=dim_list_in;
dim_list_in_max=max(dim_list_in); %max dimension to analyze
dim_list_out_max=max(dim_list_out); %max dimension to analyze
if_aug_dim=any(dim_list_out>dim_list_in);
pcon_init_method=0;
if_initpca_rot=1;
opts_pcon.if_initpca_rot=if_initpca_rot;
opts_pcon.initialize_set='pca';
opts_pcon.nshuffs=nshuffs;
opts_pcon.if_frozen=1;
opts_pcon.if_log=1;
%
%make short labels
%
dataset_labels=cell(1,nsets);
for iset=1:nsets
    dataset_labels{iset}=sets{iset}.label(1+max([find(sets{iset}.label=='/'),find(sets{iset}.label=='\')]):end);
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
ra_setup.stimulus_labels=sa_pooled.typenames;
ra_setup.sets=sets;
ra_setup.data_orig=ds;
ra_setup.metadata_orig=sas;
ra_setup.sa_consensus=sa_pooled; %metadata for ds_consensus
ra_setup.sas_components=sas_align;
%
ra_setup.nstims=nstims_all;
ra_setup.nsets=nsets;
ra_setup.nshuffs=nshuffs;
ra_setup.shuff_quantiles=shuff_quantiles; 
%
ra_setup.dim_list_in_max=dim_list_in_max;
ra_setup.dim_list_out_max=dim_list_out_max;
%
%
%do calculation for each variant of allow_scale, using same permutations
%for each shuffle
%
ra=cell(1,2);
for allow_scale=0:1
    ia=allow_scale+1;
    opts_pcon.allow_scale=allow_scale;
    %
    [ra{ia},warnings]=psg_align_stats(ds_align,sas_align,dim_list_in,dim_list_out,opts_pcon);
    %
end
%plots; display rms dev unexplained, p-values, and quantiles
ra_setup.nrows=2;
for allow_scale=0:1
    ia=allow_scale+1;
    ra_setup.row=ia;
    if (ia==1)
        figh=psg_align_stats_plot(ra{ia},ra_setup);
    else
        psg_align_stats_plot(ra{ia},setfield(ra_setup,'figh',figh));
    end
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
