% hlid_kc_knit_rs: knit several KC datasets
%
% vi: volumetric imaging (in kc_vi)
% megamat: kc_soma_nls
% tnt: controls in kc_tnt
%
%  **ask about nshuffs
%  **need to explore variants for knitted datasets from each source, and to have the datasat file names appear on figure
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
for k=1:length(sources)
    disp(sprintf('source %2.0f->%s',k,sources{k}.id));
end
source_list=getinp('choice','d',[1 length(sources)],[1:length(sources)]);
n_sources=length(source_list);
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
set_labels=cell(1,n_sources);
for source_ptr=1:n_sources
    isource=source_list(source_ptr);
    disp(sprintf(' input source type %1.0f (%s)',source_ptr,sources{isource}.id));
    opts_read.ui_filter=sources{isource}.ui_filter;
    %
    aux_read=struct;
    aux_read.nsets=0;
    aux_read.opts_check=opts_check;
    aux_read.opts_read=opts_read;
    fullnames=[];
    [d_source{source_ptr},aux_read_out]=rs_get_coordsets(fullnames,aux_read);
    nsets_each(source_ptr)=length(d_source{source_ptr}.sets);
    disp(sprintf(' %3.0f sets read.',nsets_each(source_ptr)));
    if (source_ptr==1)
        d_all=d_source{1};
    else
        d_all=rs_concat_coordsets(d_all,d_source{source_ptr});
    end
    set_labels{source_ptr}=sources{isource}.id;
end
%
%analyze with retaining stimuli that are in common to at least a criterion number (min_keep) of sources
%
if ~exist('opts_align') opts_align=struct(); end
if ~exist('opts_disp') opts_disp=struct(); end
if ~exist('opts_knit') opts_knit=struct(); end
d_align=cell(1,n_sources)
nstims_aligned=zeros(1,n_sources);
for min_keep=1:n_sources
    min_string=sprintf('min=%1.0f',min_keep);
    %
    aux_align=struct();
    aux_align.opts_check=opts_check;
    opts_align.min=min_keep;
    [d_align{min_keep},aux_align_out]=rs_align_coordsets(d_all,setfield(aux_align,'opts_align',opts_align));
    nstims_aligned(min_keep)=d_align{min_keep}.sas{1}.nstims;
    %
    if nstims_aligned(min_keep)>0
        %knit
        aux_knit=struct();
        aux_knit.opts_check=opts_check;
        opts_knit=filldefault(opts_knit,'allow_scale',1);
        opts_knit=filldefault(opts_knit,'allow_offset',1);
        opts_knit.if_pca=1;
        opts_knit.dim_max_in=5;
        opts_knit.if_dim_heuristics=1;
        opts_knit=filldefault(opts_knit,'nshuffs',50);
        opts_knit.if_stats=1;
        [d_knit,aux_knit_out]=rs_knit_coordsets(d_align{min_keep},setfield(aux_knit,'opts_knit',opts_knit));
        if opts_knit.if_stats==1
            axes('Position',[0.25,0.04,0.01,0.01]); %for text
            text(0,0,min_string,'Interpreter','none','FontSize',8);
            axis off;
            set(gcf,'Name',cat(2,get(gcf,'Name'),' ',min_string));

        end
        %plot dims 2,3,4
        for dim_select=[2:4]
            aux_disp=struct();
            opts_disp.dim_select=dim_select;
            opts_disp.coord_group_method='keeplow';
            opts_disp.data_label_setsel_method='all';
            opts_disp.connect_sets_method='all';
            opts_disp.set_labels=min_string;
            opts_disp.set_colors={'k'};
            aux_disp_out_knit=rs_disp_coordsets(d_knit,setfield(aux_disp,'opts_disp',opts_disp));
            set(gcf,'Name',cat(2,get(gcf,'Name'),' knit, ',min_string));
            %
            opts_disp_components=opts_disp;
            opts_disp_components.set_labels=set_labels;
            opts_disp_components.set_colors={'r','g','b'};
            %
            aux_disp_out_components=rs_disp_coordsets(aux_knit_out.components,setfield(aux_disp,'opts_disp',opts_disp_components));
            axes('Position',[0.25,0.04,0.01,0.01]); %for text
            text(0,0,min_string,'Interpreter','none','FontSize',8);
            axis off;
            set(gcf,'Name',cat(2,get(gcf,'Name'),' components, ',min_string));
        end
    end %nstims_aligned
end
