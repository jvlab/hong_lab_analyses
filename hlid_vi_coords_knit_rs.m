% hlid_vi_coords_knit_rs: read volumetric imaging coordinates set files,
% knit and compare
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_PCAFILT_COORDS_AUTO.
%
hlid_setup;
nsets=0;
opts_read=struct();
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
%
aux=struct;
aux.nsets=nsets;
fullnames=[];
[data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
%
aux_disp.opts_disp.set_select=1;
rs_disp_coordsets(data_read,aux_disp);
