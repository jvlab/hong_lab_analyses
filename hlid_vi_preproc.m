%hlid_vi_preproc: look at preprocessing options for KC volumetric imaging, data from George Barnum, Hong Lab
%
% 20Mar26: convert from variance ratio to F ratio
% 23Mar26: add particpation ratio for spatiotemporal part
% 31Mar26: add max_timepoints
%
%   See also:  HLID_VI_READ, HLID_VI_SPATIALFILTER, HLID_VARRATS, HLID_VI_PREPROC_PLOT, HLID_VI_VIEWPCS.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='gbarnum_mb247_soma_20241027_a_test_1.hdf5'; end
if ~exist('n_stims') n_stims=24; end
if ~exist('n_repts') n_repts=5; end
%
if ~exist('resp_measures') resp_measures={'deltaF/F','z'}; end
%
opts_read.if_remnan=1;
opts_read.if_log=0;
opts_read.if_spatialfilter=1;
fw_list=getinp('spatial kernel full width list (0 is no filter)','d',[0 Inf],[0 1 2 4 8]);
rept_list=getinp('repeat list','d',[1 n_repts],1:n_repts);
stim_list=getinp('stimulus list','d',[1 n_stims],1:n_stims);
%
opts_read.max_timepoints=getinp('max timepoints to read (0: all)','d',[0 Inf]);
%
opts_read.data_path=data_path;
opts_read.data_file=data_file;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
%
stims=hlid_vi_stimnames;
disp(data_file);
%
%
n_repts=length(rept_list);
n_stims=length(stim_list);
n_fws=length(fw_list);
n_meas=length(resp_measures);
%
var_frats=zeros(n_fws,n_meas);
part_ratios=zeros(n_fws,n_meas);
eival_sqs=zeros(n_stims*n_repts,n_fws,n_meas);
var_ratios_eacheiv=zeros(n_stims*n_repts,n_fws,n_meas); %variance ratio for each eigenvector
var_ratios_eacheiv_num=zeros(n_stims*n_repts,n_fws,n_meas); %variance ratio for each eigenvector, numerator
var_ratios_eacheiv_den=zeros(n_stims*n_repts,n_fws,n_meas); %variance ratio for each eigenvector, denominator
part_ratios_eacheiv=zeros(n_stims*n_repts,n_fws,n_meas);
part_ratios_st=zeros(n_stims*n_repts,n_fws,n_meas);
%
for fw_ptr=1:n_fws
    opts_read.sfilt_hw=fw_list(fw_ptr)/2;
    [s,opts_read_used]=hlid_vi_read(opts_read);
    n_pixels=s.n_pixels_kept;
    baseline_means=s.baseline_means;
    baseline_stdvs=s.baseline_stdvs;
    xyz=s.xyz_kept;
    if n_repts~=s.n_repts_kept
        warning('n_repts is inconsistent');
    end
    if n_stims~=s.n_stims_kept
        warning('n_stims is inconsistent');
    end
    resp_maxlength=size(s.responses,2);
    deltaF=s.responses-repmat(reshape(baseline_means,[n_pixels,1,n_repts,n_stims]),[1 resp_maxlength 1 1]);
    clear s
    %
    for meas_ptr=1:n_meas
        switch resp_measures{meas_ptr}
            case 'deltaF/F'
                v=deltaF./repmat(reshape(baseline_means,[n_pixels,1,n_repts,n_stims]),[1 resp_maxlength 1 1]);
            case 'z'
                v=deltaF./repmat(reshape(baseline_stdvs,[n_pixels,1,n_repts,n_stims]),[1 resp_maxlength 1 1]);
        end
        resp_minlength=sum(0==any(any(any(isnan(v),1),3),4));
        %
        % individual repeats
        %
        v_indiv_repts=reshape(v,[n_pixels,resp_maxlength,n_repts*n_stims]);
        v_indiv_repts=reshape(v_indiv_repts(:,[1:resp_minlength],:),[n_pixels*resp_minlength,n_repts*n_stims]);
        %
        var_rats=hlid_varrats(reshape(v_indiv_repts,[n_pixels*resp_minlength,n_repts,n_stims]));
        var_frats(fw_ptr,meas_ptr)=var_rats.frat;
        %
        %svd
        %
        [svd_u,svd_s,svd_v]=svd(v_indiv_repts,'econ'); %data=u*s*v'
        eival_sqs(:,fw_ptr,meas_ptr)=(diag(svd_s)).^2;
        part_ratios(fw_ptr,meas_ptr)=(sum(diag(svd_s)).^2)/sum(diag(svd_s).^2);
        %
        for k=1:n_stims*n_repts
            wt=svd_v(:,k); %weights for kth eigenvector (repts and stims)
            var_rats_each=hlid_varrats(reshape(wt,[1 n_repts n_stims]));
            var_ratios_eacheiv(k,fw_ptr,meas_ptr)=var_rats_each.ratio;
            var_ratios_eacheiv_num(k,fw_ptr,meas_ptr)=mean(var_rats_each.across); %for later summing
            var_ratios_eacheiv_den(k,fw_ptr,meas_ptr)=mean(var_rats_each.within); %for later summing
            %
            %svd of wts and stims to get participation ratio for repts x stims
            [svd_wu,svd_ws,svd_wv]=svd(reshape(wt,[n_repts n_stims]));
            w_dims=min(n_repts,n_stims);
            svd_ws_eivs=diag(svd_ws(1:w_dims,1:w_dims));
            part_ratios_eacheiv(k,fw_ptr,meas_ptr)=(sum(svd_ws_eivs).^2)/sum(svd_ws_eivs.^2);
            %
            %svd to get particpation ratio for spatioetemporal part
            st=svd_u(:,k); %compute particpation ratios for spatioptemporal part
            [svd_stu,svd_sts,svd_stv]=svd(reshape(st,[n_pixels resp_minlength]),'econ');
            svd_sts_eivs=diag(svd_sts(1:resp_minlength,1:resp_minlength));
            part_ratios_st(k,fw_ptr,meas_ptr)=(sum(svd_sts_eivs).^2)/sum(svd_sts_eivs.^2);
        end
    end
    clear deltaF v v_indiv_repts;
    clear svd*
    disp(sprintf(' kernel hw: %4.1f; total pixels kept: %7.0f; F ratios for dff and z: %7.4f %7.4f',opts_read.sfilt_hw,size(xyz,1),var_frats(fw_ptr,:)));
end
disp(sprintf(' suggest saving workspace as hlid_vi_preproc*.mat; can plot with hlid_vi_preproc_plot'));
