function opts_local=hlid_localopts()
% opts_local=hlid_localopts() sets local options for file names, paradigm types, etc.
%  These are typically defaults for specific file names, parameters, infixes, etc.
%
% opts_local: local options
%
% See also:  PSG_LOCALOPTS, HLID_SETUP.
%
opts_local=struct;
opts_local.coord_data_fullname_def='./data/kc_soma_nls/hlid_odor17_coords_dsid.mat'; %default full file name for a coordinate dataset
opts_local.coord_data_fullname_write_def='./data/kc_soma_nls/hlid_odor17_coords_dsid.mat'; %default full file name to write a coordinate dataset
opts_local.choice_data_fullname_def='./data/unknown_choice.mat'; %default full file name for a choice dataset
opts_local.setup_fullname_def='./data/unknown_setup.mat'; %default full file name for a setup file
%
opts_local.type_class_def='hlid'; % default type class
%
%for compatibility with built-in psg paradigms
opts_local.model_filename_def='./data/unknown_model.mat'; %default model full file name
opts_local.modeltype_def=1; %default model type parameter
opts_local.domain_list_def=cell(0);
opts_local.setup_setup_suffix=''; %suffix to convert a data file into a setup file
opts_local.domain_sigma=struct;
return
