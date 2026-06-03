%zheng_apl_read_make_consensus_rs: read, and align TNT and control files and make consensus
% no option for pca rotation, scaling, or normalization
%
%   See also: RS_GET_COORDSETS, RS_KNIT_COORDSETS, RS_GEOFIT, RS_XFORM_APPLY, RS_DISP_COORDSETS, RS_SAVEFIG.
%
%plot format customizations: to modify, set other values before running
if ~exist('dim_select_models') dim_select_models=3; end %dimensionto show in model comparison plot
if ~exist('coord_group_size_models') coord_group_size_models=3; end %coord gruop size for model_comparison plot, can be 2 or 3
if ~exist('dim_select_pa') dim_select_pa=3; end %dimensionto show in principal-axes plot
if ~exist('coord_group_size_pa') coord_group_size_pa=2; end %coord gruop size for principal-axes plot, can be 2 or 3
if ~exist('axis_view_3d') axis_view_3d=[-37.5000   30.0000]; end
if ~exist('markersize') markersize=16; end
if ~exist('linewidth') linewidth=2; end
if ~exist('callout_amount') callout_amount=0.5; end
%
if ~exist('if_sm')
    if_sm=getinp('1 to use subtract-mean datasets','d',[0 1]);
end
if if_sm
    sm_string='-sm';
    dim_reduce=1; %one less dimension
else
    sm_string=[];
    dim_reduce=0; %no dimension reduction
end
if_savemat=getinp('1 to save key variables in mat-file','d',[0 1]);
if_savefig=getinp('1 to save figure(s)','d',[0 1]);
file_name_base=cat(2,'zheng_apl_rmc',sm_string);
if (if_savemat | if_savefig)
    if_ok=0;
    while (if_ok==0)
        file_name_base=getinp('file name base','s',[],file_name_base);
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
%
%read a set to use as the adjusted dataset
data_uses.adj='TNT_label';
data_uses.ref='TNT3c';
%
data_each=struct;
%
data_use=data_uses.adj;
zheng_apl_read_align_rs;
data_each.adj=data_aligned;
clear data_aligned;
%
data_use=data_uses.ref;
zheng_apl_read_align_rs;
data_each.ref=data_aligned;
clear data_aligned;
%
%calculate consensus for adj and ref sets as
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
%clear
pcon_dim_max=nstims_all-dim_reduce;
disp(sprintf('pcon_dim_max=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',pcon_dim_max,opts_knit.pcon_init_method,opts_knit.allow_scale));
%
data_consensus=struct;
for iar=1:2
    switch iar
        case 1
            iar_string='adj';
        case 2
            iar_string='ref';
    end
    nsets=length(data_each.(iar_string).ds);
    %
    aux=struct;
    aux.opts_knit=opts_knit;
    [data_consensus.(iar_string),aux_knit_out]=rs_knit_coordsets(data_each.(iar_string),aux);
end %iar
%
%now construct geometric models
%
if ~exist('model_dim_max') model_dim_max=6; end
if (if_debug)
    if ~exist('nshuffs') nshuffs=100; end
else
    if ~exist('nshuffs') nshuffs=1000; end
end
%
model_list={'procrustes_noscale_offset'  'procrustes_scale_offset'  'affine_offset'};
colors_models={'k','b','c'};
color_tnt3c={'m'};
%
data_TNTlabel=data_consensus.adj;
data_tnt3c=data_consensus.ref;
%
aux=struct;
%
opts_geof=struct;
opts_geof.nshuffs=nshuffs;
opts_geof.model_list=model_list;
opts_geof.dimpairs_list=repmat([1:model_dim_max]',1,2);
opts_geof.dim_max_in=model_dim_max;
opts_geof.dim_max_out=model_dim_max;
opts_geof.dimpairs_method='equal'; %only analyze adj dim = ref dim
opts_geof.if_keep_transforms=1; %keep the transforms for statistical analysis
opts_geof.if_fit_summary=0;
opts_geof.if_log=0;
%
[gfs,xs,aux_out]=rs_geofit(data_TNTlabel,data_tnt3c,setfield(aux,'opts_geof',opts_geof));
%
opts_dgeo=struct;
opts_dgeo.colors_models=colors_models;
if ~exist('if_showquant') if_showquant=0; end %set to show the 0.05 quantile
opts_dgeo.if_showquant=if_showquant;
%
aux_out_dgeo=rs_disp_geofit(gfs{1}.gf,setfield(aux,'opts_dgeo',opts_dgeo));
fig_handles=aux_out_dgeo.opts_dgeo.fig_handles;
fig_names=aux_out_dgeo.opts_dgeo.fig_names;
for k=1:length(fig_handles)
    file_name_fig=cat(2,file_name_base,'_',fig_names{k});
    if (if_savefig)
        figure(fig_handles{k});
        axes('Position',[0.01,0.03,0.01,0.01]);
        text(0,0,file_name_fig,'Interpreter','none');
        axis off
        savefig(fig_handles{k},file_name_fig);
        disp(sprintf('figure saved as %s',file_name_fig));
    end
end
%
disp('Statistics')
disp(data_uses);
disp(sprintf(' nshuffs=%4.0f, if_sm=%1.0f',nshuffs,if_sm));
%
disp('Model comparison of goodness of fit via Procrustes d');
disp('    dim     d(affine)  p(affine v unif sc)   d(unif sc) p(unif sc v no sc)    d(no sc)');
disp('                        orig den  shuff den              orig den  shuff den');
ptr_affine=strmatch('affine_offset',model_list,'exact');
ptr_uniscale=strmatch('procrustes_scale_offset',model_list,'exact');
ptr_noscale=strmatch('procrustes_noscale_offset',model_list,'exact');
for k=1:model_dim_max
    d_shuff_affine_vs_uniscale=squeeze(gfs{1}.gf{k,k}.d_shuff(ptr_affine,:,ptr_uniscale,:));
    d_shuff_uniscale_vs_noscale=squeeze(gfs{1}.gf{k,k}.d_shuff(ptr_uniscale,:,ptr_noscale,:)); 
    d_affine=gfs{1}.gf{k,k}.d(ptr_affine);
    d_uniscale=gfs{1}.gf{k,k}.d(ptr_uniscale);
    d_noscale=gfs{1}.gf{k,k}.d(ptr_noscale);
    vars=[k d_affine sum(double(d_affine>=d_shuff_affine_vs_uniscale))/nshuffs d_uniscale sum(double(d_uniscale>=d_shuff_uniscale_vs_noscale))/nshuffs d_noscale];
    disp(sprintf('%8.4f   ',vars))
end
%
disp(sprintf('Axis ratio (min/max) using shuffled residuals of nested models'));
disp('                       geomean               geomean');
disp('    dim      affine    unif sc  p(unif sc)    no sc    p(no sc)');
axisrats=zeros(model_dim_max,1);
axisrats_vs_uniscale=zeros(model_dim_max,nshuffs);
axisrats_vs_noscale=zeros(model_dim_max,nshuffs);
figure; %figure for histograms of axis ratio
set(gcf,'Position',[50 100 1200 800]);
set(gcf,'Name','histograms of axis ratios');
set(gcf,'NumberTitle','off');
hist_edges=[0.3:0.02:1.0];
for k=1:model_dim_max
    A=gfs{1}.gf{k,k}.transforms{ptr_affine}.T;
    eivals=eigs(A'*A);
    axisrats(k)=sqrt(eivals(end)/eivals(1));
    for ishuff=1:nshuffs
        A_vs_uniscale=gfs{1}.gf{k,k}.transforms_nestbymodel{ptr_affine,ishuff,ptr_uniscale}.T;
        eivals_vs_uniscale=eigs(A_vs_uniscale'*A_vs_uniscale);
        axisrats_vs_uniscale(k,ishuff)=sqrt(eivals_vs_uniscale(end)/eivals_vs_uniscale(1));
        %
        A_vs_noscale=gfs{1}.gf{k,k}.transforms_nestbymodel{ptr_affine,ishuff,ptr_noscale}.T;
        eivals_vs_noscale=eigs(A_vs_noscale'*A_vs_noscale);
        axisrats_vs_noscale(k,ishuff)=sqrt(eivals_vs_noscale(end)/eivals_vs_noscale(1));
    end
    vars=[k axisrats(k) geomean(axisrats_vs_uniscale(k,:)) sum(double(axisrats(k)>=axisrats_vs_uniscale(k,:)))/nshuffs geomean(axisrats_vs_uniscale(k,:)) sum(double(axisrats(k)>=axisrats_vs_noscale(k,:)))/nshuffs];
    disp(sprintf('%8.4f   ',vars))
    if (k>1)
        subplot(2,model_dim_max-1,k-1);
        histogram(axisrats_vs_uniscale(k,:),hist_edges);
        hold on;
        plot(repmat(axisrats(k),1,2),get(gca,'YLim'),'r-');
        title(sprintf('dim %1.0f (vs unif scale)',k));
        xlabel('axis ratio')
        ylabel('counts');
        %
        subplot(2,model_dim_max-1,k-1+(model_dim_max-1));
        histogram(axisrats_vs_noscale(k,:),hist_edges);
        hold on;
        plot(repmat(axisrats(k),1,2),get(gca,'YLim'),'r-');
        title(sprintf('dim %1.0f (vs no scale)',k));
        xlabel('axis ratio')
        ylabel('counts');
    end
end
file_name_fig=cat(2,file_name_base,'_histos');
axes('Position',[0.01,0.03,0.01,0.01]);
text(0,0,file_name_fig,'Interpreter','none');
axis off
if (if_savefig)
savefig(gcf,file_name_fig);
disp(sprintf('figure saved as %s',file_name_fig));
end
%
%create transform structures from geometric models so that model
%predictions can be plotted along with control data
% 
xforms_all=struct;
model_out=struct;
for m=1:length(model_list)
    model_name=model_list{m};
    for k=1:model_dim_max
        xforms_all.(model_name).ts{1}{k}=gfs{1}.gf{k,k}.transforms{m};
    end
    xforms_all.(model_name).pipeline.desc=sprintf('transform TNT_label consensus by %s model',model_name);
    model_out.(model_name)=rs_xform_apply(data_TNTlabel,xforms_all.(model_name));
    if m==2
        data_disp=rs_concat_coordsets(model_out.(model_list{1}),model_out.(model_name));
    elseif m>=2
        data_disp=rs_concat_coordsets(data_disp,model_out.(model_name));
    end
end
data_disp=rs_concat_coordsets(data_disp,data_tnt3c);
%
% display TNT_label and models
%
opts_disp=struct();
opts_disp.set_colors=cat(2,colors_models,color_tnt3c);
opts_disp.set_labels=cat(2,model_list,'control');
opts_disp.connect_sets_linewidths=linewidth;
opts_disp.set_markersizes=markersize;
opts_disp.connect_sets_method='chain';
opts_disp.legend_interpreter='none';
opts_disp.dim_select=dim_select_models;
opts_disp.coord_group_size=coord_group_size_models;
opts_disp.callout_amount=callout_amount;
if opts_disp.coord_group_size==3
    opts_disp.axis_view=axis_view_3d;
else
    opts_disp.axis_view=2;
end
%
rs_disp_coordsets(data_disp,setfield(struct(),'opts_disp',opts_disp));
file_name_fig=cat(2,file_name_base,'_models');
axes('Position',[0.01,0.03,0.01,0.01]);
text(0,0,file_name_fig,'Interpreter','none');
axis off
if (if_savefig)
    savefig(gcf,file_name_fig);
    disp(sprintf('figure saved as %s',file_name_fig));
end
%
%create transform structures based on the principal axes
%
xforms_pa=struct;
for k=1:model_dim_max
    gfs{1}.gf{k,k}.transforms{m};
    A=gfs{1}.gf{k,k}.transforms{ptr_affine}.T;
    [v,eivals]=eigs(A'*A); %cols of v are eivecs of A'*A in matlab convention
    % A'*A*v=v*[eivals]; so (eigs in rows): v'*A*A'=[eivals]*v'
    % transpose of a column of v is a principal direction
    % mult by v on right yields the coordinates in a coordinate set of principal directions
    xforms_pa.ts{1}{k}.T=v;
    xforms_pa.ts{1}{k}.b=1;
    xforms_pa.ts{1}{k}.c=zeros(1,k);
end
model_xform=rs_xform_apply(data_disp,xforms_pa);
%
%plot in principal axes
opts_disp_pa=opts_disp;
opts_disp_pa.if_legend=-1; %extra panel with just legend
opts_disp_pa.dim_select=dim_select_pa;
opts_disp_pa.coord_group_size=coord_group_size_pa;
opts_disp_pa.axis_label_prefix='prin ax';
if opts_disp_pa.coord_group_size==3
    opts_disp_pa.axis_view=axis_view_3d;
else
    opts_disp_pa.axis_view=2;
end
rs_disp_coordsets(model_xform,setfield(struct(),'opts_disp',opts_disp_pa));
file_name_fig=cat(2,file_name_base,'_pa');
axes('Position',[0.01,0.03,0.01,0.01]);
text(0,0,file_name_fig,'Interpreter','none');
axis off
if (if_savefig)
    savefig(gcf,file_name_fig);
    disp(sprintf('figure saved as %s',file_name_fig));
end
%
%save mat file
%
if (if_savemat)
    file_name_mat=file_name_base;
    save(file_name_mat,'data_disp*','model*','gf*','opts*','-v7.3');
    disp(sprintf('variables saved as %s',file_name_mat));
end

