%hlid_setup: set up options structures for Hong Lab imaging data
%
%  Typically, run this before a psg* script, to set up opts_read and opts_plot.
%  Also sets up custom display orders for distance heatmaps.
%
% 25May24: set up custom colors for color_norays_list, so datasets can be plotted in unique colors
% 25Oct24: change plot labels to typenames, so TNT datasets will work
% 05Nov24: added display_orders.kctnt
% 06Nov24: added display_orders.kcmerge
% 14Sep25: modularize hlid_setup_func
%
% See also:  PSG_DEFOPTS, PSG_VISUALIZE_DEMO, PSG_CONSENSUS_DEMO,
% PSG_DIST_HEATMAPS_DEMO, PSG_PIPE_COORD_PROC, PSG_GEOMODELS_TEST,
% PSG_LOCALOPTS, HLID_LOCALOPTS, FILLDEFAULT.
%
h=hlid_setup_func;
hlid_opts=h.hlid_opts;
if ~exist('opts_read') opts_read=struct;end
%
if ~exist('opts_plot') opts_plot=struct; end
fns=fieldnames(h.opts_plot);
for ifn=1:length(fns)
    opts_plot=filldefault(opts_plot,fns{ifn},h.opts_plot.(fns{ifn}));
end
%
if ~exist('opts_multm_def') opts_multm_def=struct; end
fns=fieldnames(h.opts_multm_def);
for ifn=1:length(fns)
    opts_multm_def=filldefault(opts_multm_def,fns{ifn},h.opts_multm_def.(fns{ifn}));
end
disp('opts_read, opts_plot, opts_multm_def initialized for hlid data');
%
if ~exist('display_orders')
    display_orders=h.display_orders;
end
clear h
