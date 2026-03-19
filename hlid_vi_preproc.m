%hlid_vi_preproc: look at preprocessig options for KC volumetric imaging, data from George Barnum, Hong Lab
%
%   See also:  HLID_VI_READ, HLID_VI_SPATIALFILTER, HLID_VARRATS.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='gbarnum_mb247_soma_20241027_a_test_1.hdf5'; end
if ~exist('nstims') nstims=24; end
if ~exist('nrepts') nrepts=5; end
%
if ~exist('resp_measures') resp_measures={'deltaF/F','z'}; end
%
if ~exist('eigs_to_show') eigs_to_show=50; end
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
%
n_repts=length(rept_list);
n_stims=length(stim_list);
n_fws=length(fw_list);
n_meas=length(resp_measures);
%
var_ratios=zeros(n_fws,n_meas);
eival_sqs=zeros(n_stims*n_repts,n_fws,n_meas);
var_ratios_eacheiv=zeros(nstims*nrepts,n_fws,n_meas); %fariance ratio for each eigenvector
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
        var_ratios(fw_ptr,meas_ptr)=var_rats.ratio;
        %
        %svd
        %
        [svd_u,svd_s,svd_v]=svd(v_indiv_repts,'econ'); %data=u*s*v'
        eival_sqs(:,fw_ptr,meas_ptr)=(diag(svd_s)).^2;
        %
        for k=1:nstims*nrepts
            wt=svd_v(:,k); %wweights for kth eigenvector
            var_rats_each=hlid_varrats(reshape(wt,[1 n_repts n_stims]));
            var_ratios_eacheiv(k,fw_ptr,meas_ptr)=var_rats_each.ratio;
        end
    end
    clear deltaF v
    clear svd*
    disp(sprintf(' kernel hw: %4.1f; total pixels kept: %7.0f; variance ratios for dff and z: %7.4f %7.4f',opts_read.sfilt_hw,size(xyz,3),var_ratios(fw_ptr,:)));
end
%
%plot
%
if ~exist('logrange') logrange=10^2; end
eigs_to_show=min(eigs_to_show,n_repts*n_stims);
fw_labels=cell(1,n_fws);
for fw_ptr=1:n_fws
    fw_labels{fw_ptr}=sprintf('fw %2.0f',fw_list(fw_ptr));
end
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Position',[50 50 1400 800]);
set(gcf,'Name','eigenvalue anlaysis');
n_cols=3;
for meas_ptr=1:n_meas
    subplot(n_meas,n_cols,1+(meas_ptr-1)*n_cols)
    semilogy(eival_sqs(1:eigs_to_show,:,meas_ptr),'.-');
    xlabel('eigenvalue');
    ylabel(cat(2,resp_measures{meas_ptr},' var explained'));
    set(gca,'YLim',max(max(eival_sqs(:,:,meas_ptr)))*[1/logrange 1]);
    legend(fw_labels,'Location','best');
    %
    totvar=sum(eival_sqs(:,:,meas_ptr),1);
    subplot(n_meas,n_cols,2+(meas_ptr-1)*n_cols)
    semilogy(eival_sqs(1:eigs_to_show,:,meas_ptr)./repmat(totvar,eigs_to_show,1),'.-');
    xlabel('eigenvalue');
    ylabel(cat(2,resp_measures{meas_ptr},' frac var explained'));
    set(gca,'YLim',0.5*[1/logrange 1]);
    legend(fw_labels,'Location','best');
end
axes('Position',[0.01,0.01,0.01,0.01]);
text(0,0,data_file,'Interpreter','none');
axis off
