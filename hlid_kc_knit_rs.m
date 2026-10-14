% hlid_kc_knit_rs: knit several KC datasets
%
% vi: volumetric imaging (in kc_vi)
% megamat: kc_soma_nls
% tnt: controls in kc_tnt
%
%  **needs better defaults for knit stats, dims to plot
%  **need to explore variants for knitted datasets from each source
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, RS_CONCAT_COORDSETS, RS_ALIGN_COORDSETS, RS_KNIT_COORDSETS.
%
hlid_setup;
%
sources=cell(3);
sources{1}.id='megamat'; %megamat
sources{1}.ui_filter='data/kc_soma_nls/hlid_odor17_coords_kc_soma_nls_consensus*';
sources{2}.id='tnt controls';
sources{2}.ui_filter='data/kc_tnt/hlid_odor17_coords_TNT*consensus*';
sources{3}.id='volumetric';
sources{3}.ui_filter='data/kc_vi/*coords*consensus*sf0*dff*pr-all*euc_svd*';
n_sources=length(sources);
%
if ~exist('opts_read') opts_read=struct(); end
opts_read.input_type=1; %just data
opts_read.if_warn=0; %may have different sets of stimuli or stimuli in different orders
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.if_log=0;
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
%
d_source=cell(1,n_sources);
nsets_each=zeros(1,n_sources);
for isource=1:n_sources
    disp(sprintf(' input source type %1.0f (%s)',isource,sources{isource}.id));
    opts_read.ui_filter=sources{isource}.ui_filter;
    %
    aux_read=struct;
    aux_read.nsets=0;
    aux_read.opts_check=opts_check;
    aux_read.opts_read=opts_read;
    fullnames=[];
    [d_source{isource},aux_read_out]=rs_get_coordsets(fullnames,aux_read);
    nsets_each(isource)=length(d_source{isource}.sets);
    disp(sprintf(' %3.0f sets read.',nsets_each(isource)));
    if (isource==1)
        d_all=d_source{1};
    else
        d_all=rs_concat_coordsets(d_all,d_source{isource});
    end
end
aux_align=struct();
aux_align.opts_check=opts_check;
opts_align=struct();
[d_align,aux_align_out]=rs_align_coordsets(d_all,setfield(aux_align,'opts_align',opts_align));
%
aux_knit=struct();
aux_knit.opts_check=opts_check;
opts_knit=struct();
opts_knit.allow_scale=1;
opts_knit.allow_offset=1;
opts_knit.if_pca=1;
opts_knit.dim_max_in=5;
opts_knit.if_dim_heuristics=1;
opts_knit.nshuffs=50;
opts_knit.if_stats=1;
[d_knit,aux_knit_out]=rs_knit_coordsets(d_align,setfield(aux_knit,'opts_knit',opts_knit));
%
if ~exist('opts_disp') opts_disp=struct(); end
opts_disp.coord_group_method='keeplow';
opts_disp.data_label_setsel_method='all';
opts_disp.connect_sets_method='all';
aux_disp_out_knit=rs_disp_coordsets(d_knit,setfield(aux_disp,'opts_disp',opts_disp));
aux_disp_out_components=rs_disp_coordsets(aux_knit_out.components,setfield(aux_disp,'opts_disp',opts_disp));



