%zheng_apl_embed_plot: plot embeddings of TNT and control files 
%run after zheng_apl_read_align.m
%
%   See also: RS_DISP_COORDSETS.
% 
%readily-customizable plot formatting
if ~exist('axis_view') axis_view=[-37.5000   30.0000]; end
if_unif=getinp('1 for uniform coloring of each control category','d',[0 1],0);
if_indiv=getinp('1 to also plot individual datasets','d',[0 1],0);
indiv_list=getinp('list (consensus alwauys plotted)','d',[1 length(filenames.(data_use))]);
if_c2p=getinp('1 to rotate consensus into PCA space','d',[0 1],0);
%colors for plot:  first entry is for consensus
set_colors.TNT_label=opts_multm_def.color_norays_list; % {'k'  'b'  'c'  'm'  'r'  [0.7000 0.7000 0]  'g'  [0 0.5000 0]  [0.5000 0.5000 0.5000]};
set_colors.TNTin_nolabel={'k' [0.9 0 0] [0.9 0.3*(1-if_unif) 0] [0.9 0 0.5*(1-if_unif)]}; %reds
set_colors.TNT_nolabel={'k' [0.4*(1-if_unif) 0.8 0] [0 0.8 0.4*(1-if_unif)]}; %greens
set_colors.TNTin_label={'k' [0.2*(1-if_unif) 0 0.9] [0.4*(1-if_unif) 0 0.9] [0.2*(1-if_unif) 0 0.7] [0 0.2*(1-if_unif) 0.9] [0 0.4*(1-if_unif) 0.9] [0 0.2*(1-if_unif) 0.7]}; %blues
set_colors.TNT3c=cat(2,{'k'},set_colors.TNTin_nolabel(2:end),set_colors.TNT_nolabel(2:end),set_colors.TNTin_label(2:end));
set_colors.TNTall=cat(2,{'k'},set_colors.TNT_label(2:end),set_colors.TNT3c(2:end));
%
opts_disp=struct;
opts_disp.connect_sets_linewidths=2;
%
%%%%%%%%%%%%%
%
if_savemat=getinp('1 to save key variables in mat-file','d',[0 1]);
if_savefig=getinp('1 to save figure(s)','d',[0 1]);
file_name_base=cat(2,'zheng_apl_embed_plot_',data_use,sm_string);
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
if ~exist('opts_pca') opts_pca=struct(); end % for psg_pcaoffset
opts_pca=filldefault(opts_pca,'if_log',0);
opts_pca.nd_max=Inf;
%
%pcon_dim_max=getinp('maximum dimension for the consensus alignment dataset to be created','d',[1 max(dim_list_all)],pcon_dim_max);
%pcon_dim_max_comp=getinp('maximum dimension for component datasets to use (higher dimensions will be zero-padded)','d',[1 pcon_dim_max],pcon_dim_max);
%pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering','d',[-2 nsets],0);
pcon_dim_max=nstims_all-dim_reduce;
pcon_dim_max_comp=pcon_dim_max;
pcon_init_method=0;
%opts_pcon.allow_scale=getinp('1 to allow scaling for consensus','d',[0 1],opts_pcon.allow_scale);
disp(sprintf('pcon_dim_max=%3.0f, pcon_dim_max_comp=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',pcon_dim_max,pcon_dim_max_comp,pcon_init_method,opts_pcon.allow_scale));
opts_pcon.initialize_set=pcon_init_method;
opts_pcon.initialize_set='pca';
%
if if_c2p
    c2p_string='-pc';
else
    c2p_string='';
end
%
consensus=cell(pcon_dim_max,1);
z=cell(pcon_dim_max,1);
znew=cell(pcon_dim_max,1);
ts=cell(pcon_dim_max,1);
details=cell(pcon_dim_max,1);
opts_pcon_used=cell(pcon_dim_max,1);
%
ds_knitted=cell(1,pcon_dim_max); %reverse order of dimensions, 21Nov24
ds_components=cell(1,nsets); %partial datasets, aligned via Procrustes
%
for ip=1:pcon_dim_max
    z{ip}=zeros(nstims_all,ip,nsets);
    pcon_dim_use=min(ip,pcon_dim_max_comp); %pad above pcon_dim_pad
    for iset=1:nsets
        z{ip}(:,1:pcon_dim_use,iset)=ds_align{iset}{ip}(:,[1:pcon_dim_use]); %only include data up to pcon_dim_use
        z{ip}(opts_align_used.which_common_kept(:,iset)==0,:,iset)=NaN; % pad with NaN's if no data %changed from which_common to allow for more general behavior of psg_align_coordsets when opts_align.min>1
    end
    [consensus{ip},znew{ip},ts{ip},details{ip},opts_pcon_used{ip}]=procrustes_consensus(z{ip},opts_pcon);
    disp(sprintf(' creating Procrustes consensus for dim %2.0f based on datasets up to dimension %2.0f, iterations: %4.0f, final total rms dev: %8.5f',...
        ip,pcon_dim_max_comp,length(details{ip}.rms_change),sqrt(sum(details{ip}.rms_dev(:,end).^2))));
    ds_knitted{ip}=consensus{ip};
    for iset=1:nsets
        ds_components{iset}{1,ip}=znew{ip}(:,:,iset);
    end
end
%
%implement PCA rotation if requested:  note that this is applied both to consensus{ip} and to ds_components{ip}
%
ds_knitted_orig=ds_knitted;
ds_components_orig=ds_components;
if if_c2p
    for ip=1:pcon_dim_max
        knitted_centroid=mean(ds_knitted{ip},1,'omitnan');
        [ds_knitted{ip},recon_coords,var_ex,var_tot,coord_maxdiff,opts_used_pca]=psg_pcaoffset(ds_knitted{ip},knitted_centroid,opts_pca);
%        qu=opts_used_pca.qu;
%        qs=opts_used_pca.qs;
        v=opts_used_pca.qv;
        % coords=u*s*v', and recon_coords= u*s, with v'*v=I, so recon_coords=coords*v
        for iset=1:nsets
            consensus_centroid_rep=repmat(mean(ds_components{iset}{1,ip},1,'omitnan'),nstims_all,1);
            ds_components{iset}{1,ip}=consensus_centroid_rep+(ds_components{iset}{1,ip}-consensus_centroid_rep)*v(1:ip,:);
        end
    end %ip
end
%
%set up a data structure for rs_disp_coordsets, with consensus as first record
%
data_in=struct;
%
data_in.ds=cell(1,1+nsets);
data_in.sas=cell(1,1+nsets);
data_in.sets=cell(1,1+nsets);
%
data_in.ds{1}=ds_knitted;
data_in.sas{1}=setfield(sas_align{1},'metadata',[]);
sets_consensus=sets_align{1};
sets_consensus.label_long=cat(2,'consensus:',data_use,sm_string);
sets_consensus.label=sets_consensus.label;
data_in.sets{1}=sets_consensus;
%add original datasets as components
for k=1:nsets
    data_in.sets{1+k}=sets_align{k};
    data_in.sas{1+k}=sas_align{k};
    data_in.ds{1+k}=ds_components{k};
end
%opts_disp.set_alphas=0.5; only for marker-only plot
opts_disp.set_labels{1}='consensus';
for k=1:nsets
    opts_disp.set_labels{1+k}=data_in.sets{1+k}.label(strfind(data_in.sets{1+k}.label,'odor17')+7:end);
end
%
%plot consensus and individual datasets together
%
opts_disp.connect_sets_method='star';
opts_disp.axis_label_prefix=cat(2,'coord',c2p_string);
opts_disp.axis_view=axis_view;
%
if isfield(set_colors,data_use)
    opts_disp.set_colors=set_colors.(data_use);
end
opts_disp.connect_sets_color_mode='last';
aux_out=rs_disp_coordsets(data_in,setfield(struct(),'opts_disp',opts_disp));
axis_range=zeros(3,2);
axis_range(1,:)=get(gca,'XLim');
axis_range(2,:)=get(gca,'YLim');
axis_range(3,:)=get(gca,'ZLim');
file_name_fig=cat(2,file_name_base,'_all');
if (if_savefig)
    savefig(gcf,file_name_fig);
    disp(sprintf('figure saved as %s',file_name_fig));
end
set(gcf,'Name',file_name_fig);
%
%plot individual datasets
if if_indiv
    for kptr=0:length(indiv_list)
        if (kptr==0)
            k=1;
        else
            k=1+indiv_list(kptr);
        end
        opts_disp_indiv=opts_disp;
        opts_disp_indiv.axis_range='list';
        opts_disp_indiv.axis_range_list=axis_range;
        opts_disp_indiv.set_select=k;
        rs_disp_coordsets(data_in,setfield(struct(),'opts_disp',opts_disp_indiv));
        file_name_fig=cat(2,file_name_base,'_',opts_disp.set_labels{k});
        if (if_savefig)
            savefig(gcf,file_name_fig);
            disp(sprintf('figure saved as %s',file_name_fig));
        end
        set(gcf,'Name',file_name_fig);
    end
end
%
%save mat file
%
if (if_savemat)
    file_name_mat=file_name_base;
    save(file_name_mat,'data_in','opts*','-v7.3');
    disp(sprintf('variables saved as %s',file_name_mat));
end
