%zheng_apl_read_align_rs: read, and align TNT and control files
% modular version of zheng_apl_read_align
%
%  See also: RS_GET_COORDSETS.
%
if ~exist('hlid_opts')
    hlid_setup;
end
zheng_apl_files; %get file names
if ~exist('pathname') pathname='./data/kc_tnt/'; end
%
if ~exist('if_debug')
    if_debug=getinp('1 for debug mode','d',[0 1],0);
end
if if_debug
    opts_pcon.max_niters=50;
end
if ~exist('if_sm')
    if_sm=getinp('1 to use subtract-mean datasets','d',[0 1]);
end
if if_sm
    sm_string='-sm';
    dim_reduce=1; %one less dimension
else
    sm_string=[];
    dim_reduce=0; %no dimension reduction
end
%
datasets=fieldnames(filenames);
if ~exist('data_use')
    disp('data to use:')
    for k=1:length(datasets)
        disp(sprintf('%1.0f->%15s (%2.0f files)',k,datasets{k},length(filenames.(datasets{k}))));
    end
    data_use=datasets{getinp('choice','d',[1 length(datasets)])};
end
% 
nsets=length(filenames.(data_use));
fullnames=cell(1,nsets);
for k=1:nsets
    fullnames{k}=cat(2,pathname,filenames.(data_use){k},sm_string);
end
opts_read=struct();
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
%
aux=struct;
aux.nsets=nsets;
[data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
%
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
opts_align=struct();
%
[data_aligned,aux_align_out]=rs_align_coordsets(data_read,setfield(struct(),'opts_align',opts_align));
nstims_all=data_aligned.sets{1}.nstims;
disp(sprintf('total stimuli: %3.0f',nstims_all));
