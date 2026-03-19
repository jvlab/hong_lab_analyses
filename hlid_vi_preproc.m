%hlid_vi_preproc: look at preprocessig options for KC volumetric imaging, data from George Barnum, Hong Lab
%
%   See also:  HLID_VI_READ, HLID_VI_SPATIALFILTER, HLID_VARRATS.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='gbarnum_mb247_soma_20241027_a_test_1.hdf5'; end
if ~exist('nstims') nstims=24; end
if ~exist('nrepts') nrepts=5; end
%
resp_measures={'deltaF/F','z'};
%
opts_read.if_remnan=1;
opts_read.if_log=0;
opts_read.if_spatialfilter=1;
fw_list=getinp('spatial kernel full width list (0 is no filter)','d',[0 Inf],[0 1 2 4 8]);
rept_list=getinp('repeat list','d',[1 nrepts],1:nrepts);
stim_list=getinp('stimulus list','d',[1 nstims],1:nstims);
%
opts_read.data_path=data_path;
opts_read.data_file=data_file;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
%
stims=hlid_vi_stimnames;
disp(data_file);
%
var_ratios=zeros(length(fw_list),length(resp_measures));
for fw_ptr=1:length(fw_list)
    opts_read.sfilt_hw=fw_list(fw_ptr)/2;
    [s,opts_read_used]=hlid_vi_read(opts_read);
    resp_maxlength=size(s.responses,2);
    deltaF=s.responses-repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
    %
    for idffz=1:length(resp_measures)
        switch resp_measures{idffz}
            case 'deltaF/F'
                v=deltaF./repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
            case 'z'
                v=deltaF./repmat(reshape(s.baseline_stdvs,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
        end
        resp_minlength=sum(0==any(any(any(isnan(v),1),3),4));
        %
        % individual repeats
        %
        v_indiv_repts=reshape(v,[s.n_pixels_kept,resp_maxlength,s.n_repts_kept*s.n_stims_kept]);
        v_indiv_repts=reshape(v_indiv_repts(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_repts_kept*s.n_stims_kept]);
        %
        var_rats=hlid_varrats(reshape(v_indiv_repts,[s.n_pixels_kept*resp_minlength,s.n_repts_kept,s.n_stims_kept]));
        var_ratios(fw_ptr,idffz)=var_rats.ratio;
    end
    disp(sprintf(' kernel hw: %4.1f; total pixels kept: %7.0f; variance ratios for dff and z: %7.4f %7.4f',opts_read.sfilt_hw,sum(s.pixels_per_plane_kept),var_ratios(fw_ptr,:)));
end
