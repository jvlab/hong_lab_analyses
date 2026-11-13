 h% hlid_vi_coords_behav_rs: read volumetric imaging coordinates set files, compare with behavioral data from George Barnum
%
% behavioral data (valence) indicated by color
% reliability of behavior across flies indicated by size (large: less variable)
%
% 28Oct26: allow for stimuli without behavioral data
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, RS_EXTRACT_COORDSETS, HLID_VI_COORDS_KNIT_RS, HLID_VI_COORDS_KNIT_RS_AUTO, REGRESS,
%  HLID_VI_COORDS_BEHAV_AUTO_RS, HLID_VI_COORDS_BEH_DO.
%
hlid_setup;
%
if ~exist('path_beh') path_beh='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC'; end
if ~exist('filename_beh') filename_beh='attractiveness_per_fly_overview.xlsx'; end
%
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
table_beh=readtable(cat(2,path_beh,filesep,filename_beh));
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
disp(sprintf('behavioral data read from %s',filename_beh));
%
typenames_beh=table_beh.Properties.VariableDescriptions(2:end)'; %skip the index column
nstims_beh=length(typenames_beh);
rownames_beh=table_beh.fly;
rownums_beh=strmatch('202',rownames_beh); %fly names are 2024_*, etc
%
data_beh=table2array(table_beh(rownums_beh,2:end))'; %rows are stimuli, cols are preps
nstims_beh=size(data_beh,1);
npreps_beh=size(data_beh,2);
disp(sprintf('behavioral data read for %2.0f stimuli from %2.0f preps',nstims_beh,npreps_beh));
%
if ~exist('tol') tol=10^-5; end %for matching stat calcs with matlab
if ~exist('dmax') dmax=5; end
if ~exist('dims_plot') dims_plot=[2 3];end
if ~exist('line_width') line_width=2; end
if ~exist('marker_size') marker_size=12; end
if ~exist('marker_size_nobeh') marker_size_nobeh=8; end
%
dmax=getinp('maximum dimension to analyze','d',[2 10],dmax);
dims_plot=getinp('dimensions to plot (0 for none)','d',[0 dmax],dims_plot);
dims_plot=setdiff(dims_plot,[0 1]);
if ~isempty(dims_plot)
    if_plain=getinp('1 to also plot spaces without color-by-behavior','d',[0 1],0);
end
%
if ~exist('opts_read') opts_read=struct(); end
opts_read.input_type=1; %just data
opts_read.if_warn=0; %may have different sets of stimuli or stimuli in different orders
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.ui_filter='*coords*consensus*sf0*dff*pr-all*';
opts_read.if_log=0;
%
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
%
if ~exist('opts_disp') opts_disp=struct(); end
opts_disp.coord_group_method='keeplow';
opts_disp.if_legend=0;
opts_disp.callout_amount=0.5;
opts_disp.set_markersizes=marker_size;
opts_disp.set_colors='k';
%
aux=struct;
aux.nsets=0;
aux.opts_check=opts_check;
aux.opts_read=opts_read;
%
fullnames=[];
[data_read,aux_read_out]=rs_get_coordsets(fullnames,aux);
nsets=length(data_read.sets);
disp(sprintf(' %3.0f sets read.',nsets));
%
beh_string=filename_beh;
%
%calculate and plot
hlid_vi_coords_behav_do;
