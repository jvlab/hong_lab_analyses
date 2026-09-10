% hlid_vi_coords_knit_rs: read volumetric imaging coordinates set files,
% knit and compare
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_PCAFILT_COORDS_AUTO,
%   ZHENG_APL_EMBED_PLOT_RS, Rs_GET_COORDSETS, RS_KNIT_COORDSETS, RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_XFORM_SPECIFY, RS_XFORM_APPLY.
%
hlid_setup;
%
if ~exist('axis_view') axis_view=[-37.5000   30.0000]; end
if ~exist('markersize_consensus') markersize_consensus=24; end
if ~exist('markersize_component') markersize_component=16; end
if ~exist('linewidth') linewidth=2; end
%
if ~exist('prefix_remove') prefix_remove='hlid_vi_'; end
%
opts_read=struct();
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
%
aux=struct;
aux.nsets=0;
fullnames=[];
[data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
nsets=length(data_read.ds);
%
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
opts_align=struct();
%
[data_aligned,aux_align_out]=rs_align_coordsets(data_read,setfield(struct(),'opts_align',opts_align));
nstims_all=data_aligned.sets{1}.nstims;
disp(sprintf('total stimuli: %3.0f',nstims_all));
%
%aux_disp.opts_disp.set_select=1;
%rs_disp_coordsets(data_read,aux_disp);
%
if ~exist('dim_max_in') dim_max_in=10; end
if ~exist('nshuffs') nshuffs=100; end
if_stats=getinp('1 to do statistics','d',[0 1],0);
if if_stats
    dim_max_in=getinp('maximum dimension to consider','d',[1 24],dim_max_in);
    nshuffs=getinp('number of shuffles','d',[5 1000],nshuffs);
end
if_indiv=getinp('1 to also plot individual datasets','d',[0 1],0);
if if_indiv
    indiv_list=getinp('list (consensus always plotted)','d',[0 nsets]);
    indiv_list=unique([0 indiv_list]);
else
    indiv_list=[];
end
if_c2p=getinp('1 to rotate consensus into PCA space','d',[0 1],1);
%
if ~exist('dim_select_list') dim_select_list=[2:5];end
dim_select_list=getinp('dimensions to plot','d',[2 dim_max_in],dim_select_list);
%
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.if_pca=0; %if_c2p handled later
opts_knit.max_niters=1000;
opts_knit.pcon_init_method=0;
opts_knit.if_stats=if_stats;
opts_knit.dim_max_in=dim_max_in;
opts_knit.nshuffs=nshuffs;
opts_knit.if_frozen=1;
opts_knit.if_log=1;
%
opts_knit.allow_scale=getinp('1 to allow scaling','d',[0 1],opts_knit.allow_scale);
%
disp(sprintf('dim_max_in=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',opts_knit.dim_max_in,opts_knit.pcon_init_method,opts_knit.allow_scale));
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
for k=1:nsets
    data_aligned.sets{k}.label=strrep(data_aligned.sets{k}.label,prefix_remove,'');
end
[data_consensus,aux_knit_out]=rs_knit_coordsets(data_aligned,aux);
data_disp=rs_concat_coordsets(data_consensus,aux_knit_out.components); %for display: first record is components
%
opts_disp=struct;
opts_disp.connect_sets_linewidths=linewidth;
opts_disp.set_labels{1}='consensus';
%shorten the labels
for k=1:nsets
    opts_disp.set_labels{1+k}=data_read.sets{k}.label(strfind(data_read.sets{k}.label,prefix_remove)+length(prefix_remove):end);
    opts_disp.set_labels{1+k}=strrep(opts_disp.set_labels{1+k},'_','-');
end
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
    opts_disp.set_labels{1+k}=strrep(opts_disp.set_labels{1+k},'_','-');
end
if k==1
    consensus_set_label=opts_disp.set_labels{2};
else
    consensus_set_label=cat(2,'consensus (',opts_disp.set_labels{2},'...',opts_disp.set_labels{end},')');
end
%
%plot consensus and individual datasets together
%
opts_disp.set_markersizes=[markersize_consensus,repmat(markersize_component,1,nsets)];
opts_disp.connect_sets_method='star';
opts_disp.axis_label_prefix=cat(2,'coord',c2p_string);
opts_disp.axis_view=axis_view;
allow_scale_label=sprintf(' allow scaling=%1.0f',opts_knit.allow_scale);
%
%may want to control opts_disp.set_colors; first entry is for consensus
%
opts_disp.connect_sets_color_mode='last';
for dim_select=dim_select_list
    opts_disp.dim_select=dim_select;
    opts_disp.fig_name=cat(2,sprintf('dim %1.0f, ',dim_select),consensus_set_label,allow_scale_label);
    opts_disp.coord_group_method='keeplow';
    if dim_select>3
        opts_disp.if_legend=-1;
    else
        opts_disp.if_legend=1;
    end
    %
    aux_out=rs_disp_coordsets(data_disp,setfield(struct(),'opts_disp',opts_disp));
    %
    axis_range=zeros(3,2);
    axis_range(1,:)=get(gca,'XLim');
    axis_range(2,:)=get(gca,'YLim');
    axis_range(3,:)=get(gca,'ZLim');
    %
    axes('Position',[0.01,0.01,0.01,0.01]);
    text(0,0,cat(2,consensus_set_label,allow_scale_label));
    axis off
    %
    %plot individual datasets
    %
    if if_indiv
        for kptr=1:length(indiv_list)
            k=1+indiv_list(kptr);
            opts_disp_indiv=opts_disp;
            opts_disp_indiv.axis_range='list';
            opts_disp_indiv.axis_range_list=axis_range;
            opts_disp_indiv.set_select=k;
            opts_disp_indiv.fig_name=cat(2,sprintf('dim %1.0f, ',dim_select),opts_disp.set_labels{k},allow_scale_label);
            rs_disp_coordsets(data_disp,setfield(struct(),'opts_disp',opts_disp_indiv));
            %
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,cat(2,opts_disp.set_labels{k},allow_scale_label));
            axis off
        end
    end
end %plot+dim_list
