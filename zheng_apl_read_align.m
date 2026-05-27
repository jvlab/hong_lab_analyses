%zheng_apl_read_align: read, and align TNT and control files 
%
hlid_setup;
zheng_apl_files; %get file names
if ~exist('pathname') pathname='./data/kc_tnt/'; end
%
if_debug=getinp('1 for debug mode','d',[0 1],0);
if if_debug
    opts_pcon.max_niters=50;
end
if_sm=getinp('1 to use subtract-mean datasets','d',[0 1]);
if if_sm
    sm_string='-sm';
    dim_reduce=1; %one less dimension
else
    sm_string=[];
    dim_reduce=0; %no dimension reduction
end
%
datasets=fieldnames(filenames);
disp('data to use:')
for k=1:length(datasets)
    disp(sprintf('%1.0f->%15s (%2.0f files)',k,datasets{k},length(filenames.(datasets{k}))));
end
data_use=datasets{getinp('choice','d',[1 length(datasets)])};
% 
if ~exist('opts_read') opts_read=struct();end %for psg_read_coord_data
if ~exist('opts_rays') opts_rays=struct(); end %for psg_findrays
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
%
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
%
opts_align=filldefault(opts_align,'if_log',1);
%
% opts_read: options for psg_read_coorddata, can be empty or omitted
%   opts_read.if_log: 1 (default) to log (0 still shows warnings)
%   opts_read.if_warn: 1 to show warnings (defaults to 0)
%   opts_read.nfiles_max: maximum number of files to read (defaults to 100)
%   opts_read.input_type: 0 for either, 1 forces expemental data, 2 forces quadratic form, can be a scalar, or an array that is cycled through for each dataset
%   opts_read.data_fullnames: cell array of data file full names; if empty, will be requested
opts_read.data_fullnames=cell(0);
nsets=length(filenames.(data_use));
for k=1:nsets
    opts_read.data_fullnames{k}=cat(2,pathname,filenames.(data_use){k},sm_string);
end
[sets,ds,sas,rayss,opts_read_used,opts_rays_used,opts_qpred_used]=psg_get_coordsets(opts_read,opts_rays,[],nsets); %get the datasets
[sets_align,ds_align,sas_align,ovlp_array,sa_pooled,opts_align_used]=psg_align_coordsets(sets,ds,sas,opts_align); %align the stimulus names
nstims_all=sets_align{1}.nstims;
disp(sprintf('total stimuli: %3.0f',nstims_all));
