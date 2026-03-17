function [s,opts_read_used]=hlid_vi_read(opts_read)
% [s,opts_read_used]=hlid_vi_read(opts_read) reads volumetric imaging data from George Barnum, Hong Lab
% and does preliminary analysis 
%
% opts_read fields:
%  data_path: full data path
%  data_file: data file name
%  rept_list: list of repeat numbers (default, 1)
%  stim_list: list of stimulus numbers (default, 1)
%  if_log: 1 (default) to log
%  if_remnan: 1 (default) to remove pixels that are NaN in any rept or stim, -1 for just those read, 0 to keep all
%  if_keep_all_raw: 1 to keep all raw data in s, 0 (default) to delete
%
% s: structure of data read from the HDF5 file
%
% s.baseline_means(pixel_number,rept_ptr,stim_ptr): mean of raw data values during the baseline period for the response to rept_list(rept_ptr) and stim_list(stim_ptr)
% s.baseline_stdvs(pixel_number,rept_ptr,stim_ptr): standard dev of raw data values during the baseline period for the response to rept_list(rept_ptr) and stim_list(stim_ptr)
% s.responses(pixel_number,frame,rept_ptr,stim_ptr): raw data values for the response to rept_list(rept_ptr) and stim_list(stim_ptr)
%
% opts_read_used: options used
%
%   See also:  HLID_VI_EXPLORE, H5INFO, H5READ.
%
if nargin<1
    opts_read=struct;
end
opts_read=filldefault(opts_read,'data_path','C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\');
opts_read=filldefault(opts_read,'data_file','gbarnum_mb247_soma_20241027_a_test.hdf5');
opts_read=filldefault(opts_read,'rept_list',1);
opts_read=filldefault(opts_read,'stim_list',1);
opts_read=filldefault(opts_read,'if_log',1);
opts_read=filldefault(opts_read,'if_remnan',1);
opts_read=filldefault(opts_read,'if_keep_all_raw',0);
opts_read_used=opts_read;
%
s=struct;
s.opts_read=opts_read;
fullname=cat(2,opts_read.data_path,filesep,opts_read.data_file);
s.info=h5info(fullname);
s.n_planes=double(h5read(fullname,'/n_planes'));
s.odor_ids=double(h5read(fullname,'/odor_ids'));
s.xyz_all=[int16(h5read(fullname,'/x')),int16(h5read(fullname,'/y')),int16(h5read(fullname,'/z'))]; %use int16 to save space
%
s.n_pixels_all=size(s.xyz_all,1);
s.onset_indexes_all=h5read(fullname,'/onset_indexes');
s.offset_indexes_all=h5read(fullname,'/offset_indexes');
s.n_repts_all=size(s.onset_indexes_all,1);
s.n_stims_all=size(s.offset_indexes_all,2);
s.n_repts_kept=length(opts_read.rept_list);
s.n_stims_kept=length(opts_read.stim_list);
s.onset_indexes_kept=s.onset_indexes_all(opts_read.rept_list,opts_read.stim_list);
s.offset_indexes_kept=s.offset_indexes_all(opts_read.rept_list,opts_read.stim_list);
%
info_pv=h5info(fullname,'/pixel_values');
s.n_timepoints=info_pv.Dataspace.MaxSize(2);
n_pixels_check=info_pv.Dataspace.MaxSize(1);
if n_pixels_check~=s.n_pixels_all
    warning(sprintf('mismatch between expected %10.0f and found number %10.0f number of pixels',s.n_pixels_all,n_pixels_check));
end
%read pixel values
%
if opts_read.if_remnan==1
    pixel_values=h5read(fullname,'/pixel_values'); %read all repts and stims
    nan_pixels=any(any(any(isnan(pixel_values),2),3),4);
    if opts_read.if_log
        disp(sprintf('%10.0f pixel values in %7.0f pixels across all repts and stims; %7.0f have a NaN in some frame',...
            numel(pixel_values),size(pixel_values,1),sum(nan_pixels)));
    end
    pixels_keep=find(nan_pixels==0);
    pixel_values=pixel_values(:,:,opts_read.rept_list,opts_read.stim_list); %pixel valuse in requested repts and stims
else
    pixel_values=zeros(s.n_pixels_all,s.n_timepoints,s.n_repts_kept,s.n_stims_kept); %only read the repts and stims rquested
    for rept_ptr=1:s.n_repts_kept
        for stim_ptr=1:s.n_stims_kept
            rept=opts_read.rept_list(rept_ptr);
            stim=opts_read.stim_list(stim_ptr);
            pixel_values(:,:,rept_ptr,stim_ptr)=h5read(fullname,'/pixel_values',[1 1 rept stim],[s.n_pixels_all s.n_timepoints 1 1]);
        end
    end
    nan_pixels=any(any(any(isnan(pixel_values),2),3),4);
    if opts_read.if_log
        disp(sprintf('%10.0f pixel values in %7.0f pixels across specified repts and stims; %7.0f have a NaN in some frame',...
            numel(pixel_values),size(pixel_values,1),sum(nan_pixels)));
    end
    if opts_read.if_remnan==-1
        pixels_keep=find(nan_pixels==0);
    else
        pixels_keep=[1:size(pixel_values,1)];
    end
end
s.n_pixels_kept=length(pixels_keep);
s.xyz_kept=s.xyz_all(pixels_keep,:);
s.plane_list_all=double(unique(s.xyz_all(:,3)));
s.n_planes_with_data_all=length(s.plane_list_all);
s.pixels_per_plane_all=zeros(1,s.n_planes_with_data_all);
s.pixels_per_plane_all_kept=zeros(1,s.n_planes_with_data_all); %the count of the no-Nan pixels in all planes prior to removal of NaN, but may include planes in wich all pixels have been removed
for plane_ptr=1:s.n_planes_with_data_all
    plane=s.plane_list_all(plane_ptr);
    s.pixels_per_plane_all(plane_ptr)=sum(sum(s.xyz_all(:,3)==plane));
    s.pixels_per_plane_all_kept(plane_ptr)=sum(sum(s.xyz_all(pixels_keep,3)==plane));
end
s.pixel_values_kept=pixel_values(pixels_keep,:,:,:);
if opts_read.if_log
    for plane_ptr=1:s.n_planes_with_data_all
        plane=s.plane_list_all(plane_ptr);
        disp(sprintf(' plane %3.0f: pixels present: %7.0f, pixels kept: %7.0f',...
            s.plane_list_all(plane_ptr),s.pixels_per_plane_all(plane_ptr),s.pixels_per_plane_all_kept(plane_ptr)));
    end
    disp(sprintf('      total pixels present: %7.0f, pixels kept: %7.0f',sum(s.pixels_per_plane_all),sum(s.pixels_per_plane_all_kept)));
    disp(sprintf(' ck:  total pixels present: %7.0f, pixels kept: %7.0f',size(pixel_values,1),length(pixels_keep)))
end
s.plane_list_kept=s.plane_list_all(s.pixels_per_plane_all_kept>0);
s.n_planes_with_data_kept=length(s.plane_list_kept);
s.pixels_per_plane_kept=s.pixels_per_plane_all_kept(s.pixels_per_plane_all_kept>0);
%
%Identify and extract baseline and stim values
%
% From George (with note that the last value of the range is to be excluced, Python notation)
%To get the volume level intervals, the baseline period should be [1:floor(onset_index/n_planes)+1] and the stimulus period should be [ceil(onset_index/n_planes)+1:floor(offset_index/n_planes)+1].
%This should give all the volumes entirely contained within the relevant period, and exclude the volume that straddles the onset of the stimulus.
%
%
s.baseline_frame_range=cat(3,ones(s.n_repts_kept,s.n_stims_kept),floor(s.onset_indexes_kept/s.n_planes));
s.response_frame_range=cat(3,ceil(s.onset_indexes_kept/s.n_planes)+1,floor(s.offset_indexes_kept./s.n_planes));
s.response_lengths=1+diff(s.response_frame_range,1,3);
response_max_length=max(s.response_lengths(:));
%
s.baseline_means=zeros(s.n_pixels_kept,s.n_repts_kept,s.n_stims_kept);
s.baseline_stdvs=zeros(s.n_pixels_kept,s.n_repts_kept,s.n_stims_kept);
s.responses=NaN(s.n_pixels_kept,response_max_length,s.n_repts_kept,s.n_stims_kept);
%
for rept_ptr=1:s.n_repts_kept
    for stim_ptr=1:s.n_stims_kept
        baseline_frame_lims=squeeze(s.baseline_frame_range(rept_ptr,stim_ptr,:));
        baseline_frames=baseline_frame_lims(1):baseline_frame_lims(2);
        s.baseline_means(:,rept_ptr,stim_ptr)=mean(s.pixel_values_kept(:,baseline_frames,rept_ptr,stim_ptr),2,'omitnan');
        s.baseline_stdvs(:,rept_ptr,stim_ptr)=std(s.pixel_values_kept(:,baseline_frames,rept_ptr,stim_ptr),0,2,'omitnan');
        response_frames=[s.response_frame_range(rept_ptr,stim_ptr,1):s.response_frame_range(rept_ptr,stim_ptr,2)];
        s.responses(:,[1:s.response_lengths(rept_ptr,stim_ptr)],rept_ptr,stim_ptr)=s.pixel_values_kept(:,response_frames,rept_ptr,stim_ptr);
    end
end
if opts_read.if_keep_all_raw==0
    s=rmfield(s,'pixel_values_kept');
end
return
end

