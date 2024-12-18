%hlid_read_coorddata_demo: read a calcium imaging data coordinate file suitable for the psg package
%
%  See also:  HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO, PSG_GET_COORDSETS, PSG_READ_COORDDATA.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
if ~exist('opts_read') opts_read=struct;end
opts_read=filldefault(opts_read,'if_log',1);
opts_read=filldefault(opts_read,'if_justsetup',0);
opts_read=filldefault(opts_read,'data_fullname_def',hlid_opts.coord_data_fullname_def);
opts_read=filldefault(opts_read,'setup_fullname_def',[]);
opts_read=filldefault(opts_read,'domain_list',hlid_opts.domain_list_def);
opts_read=filldefault(opts_read,'setup_suffix',hlid_opts.setup_setup_suffix);
opts_read=filldefault(opts_read,'type_class_aux',hlid_opts.type_class_def);
%
[d,sa,opts_read_used,pipeline]=psg_read_coorddata([],[],opts_read);
disp('read a single dataset');
nsets=getinp('number of sets to read together','d',[1 10],2);
[sets,ds,sas,rayss,opts_read_used2,opts_rays_used,opts_qpred_used]=psg_get_coordsets(setfield(opts_read,'input_type',1),[],[],nsets);
disp('read multiple datasets')
