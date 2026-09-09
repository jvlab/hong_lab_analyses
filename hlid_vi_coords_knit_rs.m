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
if_indiv=getinp('1 to also plot individual datasets','d',[0 1],0);
if if_indiv
    indiv_list=getinp('list (consensus always plotted)','d',[1 nsets]);
else
    indiv_list=[];
end
if_c2p=getinp('1 to rotate consensus into PCA space','d',[0 1],0);
%colors for plot:  first entry is for consensus
%set_colors=[];
%
opts_disp=struct;
opts_disp.connect_sets_linewidths=linewidth;
%
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.max_niters=1000;
opts_knit.pcon_init_method=0;
opts_knit.if_stats=0; %no statistics
opts_knit.if_frozen=1;
opts_knit.if_log=1;
%
opts_knit.allow_scale=getinp('1 to allow scaling','d',[0 1],opts_knit.allow_scale);
%
if ~exist('pcon_dim_max') pcon_dim_max=10; end
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
    opts_disp.set_labels{1+k}=data_disp.sets{1+k}.label(strfind(data_disp.sets{1+k}.label,'hlid_vi_')+8:end);
    opts_disp.set_labels{1+k}=strrep(opts_disp.set_labels{1+k},'_','-');
end
%
%plot consensus and individual datasets together
%
opts_disp.set_markersizes=[markersize_consensus,repmat(markersize_component,1,nsets)];
opts_disp.connect_sets_method='star';
opts_disp.axis_label_prefix=cat(2,'coord',c2p_string);
opts_disp.axis_view=axis_view;
%
%if isfield(set_colors,data_use)
%    opts_disp.set_colors=set_colors.(data_use);
%end
opts_disp.connect_sets_color_mode='last';
aux_out=rs_disp_coordsets(data_disp,setfield(struct(),'opts_disp',opts_disp));
axis_range=zeros(3,2);
axis_range(1,:)=get(gca,'XLim');
axis_range(2,:)=get(gca,'YLim');
axis_range(3,:)=get(gca,'ZLim');
%
%plot individual datasets
%
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
        % file_name_fig=cat(2,file_name_base,'_',opts_disp.set_labels{k});
        % if (if_savefig)
        %     savefig(gcf,file_name_fig);
        %     disp(sprintf('figure saved as %s',file_name_fig));
        % end
        % set(gcf,'Name',file_name_fig);
    end
end
