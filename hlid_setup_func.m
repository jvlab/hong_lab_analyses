function h=hlid_setup_func
%h=hlid_setup_func returns default optoins for hlid datasets
%
%  Typically, run hlid_setup, which calls this, before a psg* script, to set up opts_read and opts_plot.
%  Also sets up custom display orders for distance heatmaps.
%
% See also:  HLID_SETUP, HLID_LOCALOPTS.
%
h=struct;
%
h.hlid_opts=hlid_localopts();
%
h.opts_read=struct;
h.opts_read.if_log=1;
h.opts_read.if_justsetup=0;
h.opts_read.data_fullname_def=h.hlid_opts.coord_data_fullname_def;
h.opts_read.setup_fullname_def=[];
h.opts_read.domain_list=h.hlid_opts.domain_list_def;
h.opts_read.setup_suffix=h.hlid_opts.setup_setup_suffix;
h.opts_read.type_class_aux=h.hlid_opts.type_class_def;
h.opts_read.if_spray=1; %rays can be defined by single points
h.opts_read.if_data_only=1; %no "model" type for input data
%
h.opts_plot=struct;
h.opts_plot.if_use_rays=0; %ray structure is not used for plotting
h.opts_plot.label_sets=1; %label the first set when plotting
h.opts_plot.label_list='typenames';
%opts_plot=filldefault(opts_plot,'label_list',...
%    {'1-5ol','1-6ol','1-8ol','2-but','2h','6al','B-cit','IaA','Lin','aa','benz','eb','ep','ms','pa','t2h','va'});
%
h.opts_multm_def=struct;
h.opts_multm_def.color_norays_list={'k','b','c','m','r',[0.7 0.7 0],'g',[0 0.5 0],[0.5 0.5 0.5]}; %colors for individual datasets
%
h.display_orders=struct;
h.display_orders.kcclust={'2h','IaA','pa','2-but','eb','ep','aa','va','B-cit','Lin','6al','t2h','1-8ol','1-5ol','1-6ol','benz','ms'};
h.display_orders.kctnt={'1-5ol','1-6ol','1-8ol','2-but','2h','6al','IaA','PAA','eb','ep','ms','pa','t2h','B-myr','euc','-aPine','pfo'};
h.display_orders.kcmerge={'2h','IaA','pa','2-but','eb','ep','6al','t2h','1-8ol','1-5ol','1-6ol','PAA','ms','B-myr','euc','-aPine','pfo'};
return
end
