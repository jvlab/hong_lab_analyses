% hlid_vi_coords_knit_auto_rs: read volumetric imaging coordinates set files,
% automated knit and compare, based on outputs of hlid_vi_pcafilt_coords_auto.
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_PCAFILT_COORDS_AUTO,
%   HLID_VI_COORDS_KNIT_RS,
%   ZHENG_APL_EMBED_PLOT_RS, Rs_GET_COORDSETS, RS_KNIT_COORDSETS, RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_XFORM_SPECIFY, RS_XFORM_APPLY.
%
hlid_setup;
%
if ~exist('axis_view') axis_view=[-37.5000   30.0000]; end
if ~exist('markersize_consensus') markersize_consensus=24; end
if ~exist('markersize_component') markersize_component=16; end
if ~exist('linewidth') linewidth=2; end
if ~exist('callout_amount') callout_amount=0.5; end
%
if ~exist('prefix_remove') prefix_remove='hlid_vi_'; end
%
opts_read=struct();
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.if_log=0;
%
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
opts_align.if_log=0;
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
%
if ~exist('results_file') results_file='hlid_vi_pcafilt_coords_auto_11Sep26.mat'; end
if ~exist('coord_path') coord_path='./data/kc_vi'; end
results_file=getinp('file with results from hlid_vi_pcafilt_coords_auto','s',[],results_file);
coord_path=getinp('path to coordinate files created by hlid_vi_pcafilt_coords_auto','s',[],coord_path);
cord_path=strrep(strrep(coord_path,'/',filesep),'\',filesep);
%
load(results_file,'results');
results_dims=size(results);
nsets=results_dims(1);
nvariants=prod(results_dims(2:end));
results_res=reshape(results,nsets,nvariants);
disp('loaded results, with dimensions')
disp(results_dims)
disp(sprintf('nsets=%3.0f  nvariants=%5.0f',nsets,nvariants));
disp(sprintf('first coord file of first variant: %s',results_res{1,1}.coord_file));
disp(sprintf(' last coord file of first variant: %s',results_res{end,1}.coord_file));
disp(sprintf('first coord file of  last variant: %s',results_res{1,end}.coord_file));
disp(sprintf(' last coord file of  last variant: %s',results_res{end,end}.coord_file));
%
if ~exist('dim_max_in') dim_max_in=10; end
if ~exist('nshuffs') nshuffs=100; end
%
dim_max_in=getinp('maximum dimension to consider','d',[1 24],dim_max_in);
nshuffs=getinp('number of shuffles','d',[0 1000],nshuffs);
%
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.if_pca=0; %if_c2p handled later
opts_knit.max_niters=1000;
opts_knit.pcon_init_method=0;
opts_knit.if_stats=1;
opts_knit.if_plot=0;
opts_knit.dim_max_in=dim_max_in;
opts_knit.nshuffs=nshuffs;
opts_knit.if_frozen=1;
opts_knit.if_log=0;
%
results_knit=cell(1,nvariants);
for ivariant=1:nvariants
    fullnames=cell(1,nsets);
    for iset=1:nsets
        fullnames{iset}=cat(2,coord_path,filesep,results_res{iset,ivariant}.coord_file);
    end
    aux=struct;
    aux.nsets=nsets;
    [data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
    disp('***********')
    disp(sprintf('variant %5.0f',ivariant));
    disp(sprintf('first file %s',fullnames{1}));
    disp(sprintf(' last file %s',fullnames{nsets}));
    %
    [data_aligned,aux_align_out]=rs_align_coordsets(data_read,setfields(struct(),{'opts_align','opts_check'},{opts_align,opts_check}));
    nstims_all=data_aligned.sets{1}.nstims;
    disp('aligned');
    %
    aux.opts_knit=opts_knit;
    [data_consensus,aux_knit_out]=rs_knit_coordsets(data_aligned,setfields(struct(),{'opts_knit','opts_check'},{opts_knit,opts_check}));
    disp(sprintf('knit with dim_max_in=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',opts_knit.dim_max_in,opts_knit.pcon_init_method,opts_knit.allow_scale));
    %collect key rmsw variance and shuffle stats
    results_knit{ivariant}.rmsavail_overall=aux_knit_out.knit_stats.rmsavail_overall;
    results_knit{ivariant}.rmsdev_overall=aux_knit_out.knit_stats.rmsdev_overall;
    results_knit{ivariant}.rmsdev_overall_shuff=reshape(aux_knit_out.knit_stats.rmsdev_overall_shuff(:,1,1,:,1),[dim_max_in nshuffs]); %d1: dimension, d2: which shuffle
end
