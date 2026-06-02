%zheng_apl_embed_plot: plot embeddings of TNT and control files 
%run after zheng_apl_read_align_rs.m
%
%   See also: RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_XFORM_SPECIFY, RS_XFORM_APPLY.
% 
%readily-customizable plot formatting; these variables can be set before running to override the defaults
if ~exist('axis_view') axis_view=[-37.5000   30.0000]; end
if ~exist('markersize_consensus') markersize_consensus=24; end
if ~exist('markersize_component') markersize_component=16; end
if ~exist('linewidth') linewidth=2; end
%
if_unif=getinp('1 for uniform coloring of each control category','d',[0 1],0);
if_indiv=getinp('1 to also plot individual datasets','d',[0 1],0);
if if_indiv
    indiv_list=getinp('list (consensus always plotted)','d',[1 length(filenames.(data_use))]);
else
    indiv_list=[];
end
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
opts_disp.connect_sets_linewidths=linewidth;
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
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.max_iters=1000; %nonstandard max
opts_knit.pcon_init_method=0;
opts_knit.if_stats=0; %no statistics
opts_knit.if_frozen=1;
opts_knit.if_log=1;
%
pcon_dim_max=nstims_all-dim_reduce;
disp(sprintf('pcon_dim_max=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',pcon_dim_max,opts_knit.pcon_init_method,opts_knit.allow_scale));
%
if if_c2p
    c2p_string='-pc';
else
    c2p_string='';
end
%
%knit to find consensus
%
aux=struct;
aux.opts_knit=opts_knit;
[data_consensus,aux_knit_out]=rs_knit_coordsets(data_aligned,aux);
data_disp=rs_concat_coordsets(data_consensus,aux_knit_out.components); %for display: first record is components
%
%PCA rotation about centroid if requested
%
if if_c2p
    opts_xform=struct();
    opts_xform.mode='offset_pca';
    opts_xform.source='local';
    opts_xform.centering_specifier='centroid';
    xforms=rs_xform_specify(data_consensus,setfield(struct(),'opts_xform',opts_xform));
    data_consensus_orig=data_consensus;
    data_consensus=rs_xform_apply(data_consensus_orig,xforms);
    data_disp_orig=data_disp;
    data_disp=rs_xform_apply(data_disp_orig,xforms);
end
%
%display options
%
opts_disp.set_labels{1}='consensus';
for k=1:nsets
    opts_disp.set_labels{1+k}=data_disp.sets{1+k}.label(strfind(data_disp.sets{1+k}.label,'odor17')+7:end);
end
%
%plot consensus and individual datasets together
%
opts_disp.set_markersizes=[markersize_consensus,repmat(markersize_component,1,nsets)];
opts_disp.connect_sets_method='star';
opts_disp.axis_label_prefix=cat(2,'coord',c2p_string);
opts_disp.axis_view=axis_view;
%
if isfield(set_colors,data_use)
    opts_disp.set_colors=set_colors.(data_use);
end
opts_disp.connect_sets_color_mode='last';
aux_out=rs_disp_coordsets(data_disp,setfield(struct(),'opts_disp',opts_disp));
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
        rs_disp_coordsets(data_disp,setfield(struct(),'opts_disp',opts_disp_indiv));
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
    save(file_name_mat,'data_disp*','opts*','-v7.3');
    disp(sprintf('variables saved as %s',file_name_mat));
end
