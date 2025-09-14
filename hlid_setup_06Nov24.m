%hlid_setup: set up options structures for Hong Lab imaging data
%
%  Typically, run this before a psg* script, to set up opts_read and opts_plot.
%  Also sets up custom display orders for distance heatmaps.
%
% 25May24: set up custom colors for color_norays_list, so datasets can be plotted in unique colors
% 25Oct24: change plot labels to typenames, so TNT datasets will work
% 05Nov24: added display_orders.kctnt
% 06Nov24: added display_orders.kcmerge
%
% See also:  PSG_DEFOPTS, PSG_VISUALIZE_DEMO, PSG_CONSENSUS_DEMO,
% PSG_DIST_HEATMAPS_DEMO, PSG_PIPE_COORD_PROC, PSG_GEOMODELS_TEST,
% PSG_LOCALOPTS, HLID_LOCALOPTS, FILLDEFAULT.
%
hlid_opts=hlid_localopts;
if ~exist('opts_read') opts_read=struct;end
opts_read=filldefault(opts_read,'if_log',1);
opts_read=filldefault(opts_read,'if_justsetup',0);
opts_read=filldefault(opts_read,'data_fullname_def',hlid_opts.coord_data_fullname_def);
opts_read=filldefault(opts_read,'setup_fullname_def',[]);
opts_read=filldefault(opts_read,'domain_list',hlid_opts.domain_list_def);
opts_read=filldefault(opts_read,'setup_suffix',hlid_opts.setup_setup_suffix);
opts_read=filldefault(opts_read,'type_class_aux',hlid_opts.type_class_def);
opts_read=filldefault(opts_read,'if_spray',1); %rays can be defined by single points
opts_read=filldefault(opts_read,'if_data_only',1); %no "model" type for input data
%
if ~exist('opts_plot') opts_plot=struct; end
opts_plot=filldefault(opts_plot,'if_use_rays',0); %ray structure is not used for plotting
opts_plot=filldefault(opts_plot,'label_sets',1); %label the first set when plotting
opts_plot=filldefault(opts_plot,'label_list','typenames');
%opts_plot=filldefault(opts_plot,'label_list',...
%    {'1-5ol','1-6ol','1-8ol','2-but','2h','6al','B-cit','IaA','Lin','aa','benz','eb','ep','ms','pa','t2h','va'});
%
if ~exist('opts_multm_def') opts_multm_def=struct; end
opts_multm_def=filldefault(opts_multm_def,'color_norays_list',{'k','b','c','m','r',[0.7 0.7 0],'g',[0 0.5 0],[0.5 0.5 0.5]}); %colors for individual datasets
disp('opts_read, opts_plot, opts_multim_def initialized for hlid data');
%
if ~exist('display_orders')
    display_orders=struct;
    display_orders.kcclust={'2h','IaA','pa','2-but','eb','ep','aa','va','B-cit','Lin','6al','t2h','1-8ol','1-5ol','1-6ol','benz','ms'};
    display_orders.kctnt={'1-5ol','1-6ol','1-8ol','2-but','2h','6al','IaA','PAA','eb','ep','ms','pa','t2h','B-myr','euc','-aPine','pfo'};
    display_orders.kcmerge={'2h','IaA','pa','2-but','eb','ep','6al','t2h','1-8ol','1-5ol','1-6ol','PAA','ms','B-myr','euc','-aPine','pfo'};
end
