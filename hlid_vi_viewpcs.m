%hlid_vi_viewpcs: view principal components of KC volumetric imaging, data from George Barnum, Hong Lab
%
%   See also:  HLID_VI_READ, HLID_VI_SPATIALFILTER, HLID_VARRATS, HLID_VI_PREPROC, HLID_VI_VIEWPCS_UTIL.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='gbarnum_mb247_soma_20241027_a_test_1.hdf5'; end
if ~exist('n_stims') n_stims=24; end
if ~exist('n_repts') n_repts=5; end
if ~exist('n_tbins') n_tbins=5; end %default number of time bins to use
if ~exist('n_stpcs') n_stpcs=2; end % number of pc's to show of each spatiotemporal PC
%
if ~exist('resp_measures') resp_measures={'deltaF/F','z'}; end
%
opts_read.if_remnan=1;
opts_read.if_log=0;
for k=1:length(resp_measures)
    disp(sprintf('%1.0f->response measure %s',k,resp_measures{k}));
end
resp_measure=resp_measures{getinp('choice','d',[1 length(resp_measures)],1)};
opts_read.if_spatialfilter=getinp('1 for spatial filter','d',[0 1],0);
if opts_read.if_spatialfilter
    sfilt_fw=getinp('spatial kernel full width (0 is no filter)','d',[0 Inf],1);
    opts_read.sfilt_hw=sfilt_fw/2;
    sf_string=sprintf('sf: kernel hw=%3.1f',opts_read.sfilt_hw);
else
    sf_string='sf: none';
end
rept_list=getinp('repeat list','d',[1 n_repts],1:n_repts);
stim_list=getinp('stimulus list','d',[1 n_stims],1:n_stims);
n_repts=length(rept_list);
n_stims=length(stim_list);
%
stims=hlid_vi_stimnames;
disp(sprintf('reading %s',data_file));
read_data_file_short=strrep(data_file,'.hdf5','');
rept_string=cat(2,'repts: ',sprintf(' %2.0f',rept_list));
%
opts_read.data_path=data_path;
opts_read.data_file=data_file;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
[s,opts_read_used]=hlid_vi_read(opts_read);
%
tstring=cat(2,read_data_file_short,':',sf_string,',',resp_measure,', ',rept_string);
disp(sprintf('read %s',tstring));
%
n_pixels=s.n_pixels_kept;
baseline_means=s.baseline_means;
baseline_stdvs=s.baseline_stdvs;
xyz=s.xyz_kept;
planes=unique(xyz(:,3));
n_planes=length(planes);
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
switch resp_measure
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
%svd
%
disp('computing pcs')
[svd_u,svd_s,svd_v]=svd(v_indiv_repts,'econ'); %data=u*s*v'
%
clear v v_indiv_repts deltaF baseline_means baseline_stdvs
%
if_done=0;
ipc=[];
if_hv=1;
if_unifscale=1;
svd_v_max=max(abs(svd_v(:)));
svd_u_max=max(abs(svd_u(:)));
svd_s_dsq=diag(svd_s).^2;
while (if_done==0)
    ipc=getinp('pc to show (0 to end)','d',[0 n_repts*n_stims],ipc);
    if ipc==0
        if_done=1;
    else
        pc_string=sprintf('pc %3.0f',ipc);
        pc_string2=cat(2,pc_string,sprintf(' frac of variance: %8.6f, rel. frac vs. pc1: %8.6f',...
            svd_s_dsq(ipc)/sum(svd_s_dsq),svd_s_dsq(ipc)/svd_s_dsq(1)));
        n_tbins=getinp('number of time bins','d',[1 resp_minlength],n_tbins);
        tbins=repmat(floor(resp_minlength/n_tbins),1,n_tbins);
        if_ok=0;
        while (if_ok==0)
            tbins=getinp(sprintf('time bins (sums to <=%2.0f)',resp_minlength),'d',[1 resp_minlength],tbins);
            if_ok=double(sum(tbins)<=resp_minlength) & (length(tbins)==n_tbins);
        end
        if_hv=getinp('1 for horizontal, 2 for vertical','d',[1 2],if_hv);
        if_unifscale=getinp('1 for uniform scale, 0 to scale each to max','d',[0 1],if_unifscale);
        opts_viewpcs=struct;
        opts_viewpcs.if_hv=if_hv;
        opts_viewpcs.if_unifscale=if_unifscale;
        %
        n_stpcs=getinp('number of spatiotemporal pcs to show','d',[0 n_stpcs],n_stpcs);
        %
        bin_ranges=[(1+[0 cumsum(tbins(1:end-1))]);cumsum(tbins)];
        spatem=reshape(svd_u(:,ipc),[n_pixels resp_minlength]);
        repstm=reshape(svd_v(:,ipc),[n_repts n_stims]);
        spatem_binned=zeros(n_pixels,n_tbins);
        for tbin=1:n_tbins
            spatem_binned(:,tbin)=mean(spatem(:,[bin_ranges(1,tbin):bin_ranges(2,tbin)]),2);
        end
        %
        %plot spatiotemporal component with a plane in each row
        %
        figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,pc_string,' spatiotemp detail, ',tstring));
        set(gcf,'Position',[50 50 1400 800]);
        %
        opts_viewpcs.bin_ranges=bin_ranges;
        opts_viewpcs.n_cols=n_tbins;
        opts_viewpcs.col_off=0;
        hlid_vi_viewpcs_util(xyz,spatem_binned,opts_viewpcs);
        %
        axes('Position',[0.01,0.04,0.01,0.01]);
        text(0,0,pc_string2,'Interpreter','none');
        axis off
        axes('Position',[0.01,0.01,0.01,0.01]);
        text(0,0,tstring,'Interpreter','none');
        axis off
        %
        %plot summary of spatiotemporal component 
        %
        figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,pc_string,' summary, ',tstring));
        set(gcf,'Position',[50 50 1400 800]);
        spatem_avg=mean(spatem,2);
        %
        n_sumcols=2*(1+n_stpcs);
        opts_viewpcs.bin_ranges=[1 resp_minlength]';
        opts_viewpcs.n_cols=n_sumcols;
        opts_viewpcs.col_off=0;
        handles=hlid_vi_viewpcs_util(xyz,spatem_avg,opts_viewpcs);
        axpos=handles{1}.Position;
        axes('Position',[axpos(1) axpos(2)+axpos(4)+0.02,0.01,0.01]);
        text(0,0,'mean');
        axis off
        if n_stpcs>0
            %svd of spatiotemporal part
            [svd_stu,svd_sts,svd_stv]=svd(spatem,'econ');
            svd_sts_dsq=diag(svd_sts).^2;
            for istpc=1:n_stpcs
                handles=hlid_vi_viewpcs_util(xyz,svd_stu(:,istpc),setfield(opts_viewpcs,'col_off',istpc));
                axpos=handles{1}.Position;
                axes('Position',[axpos(1) axpos(2)+axpos(4)+0.02,0.01,0.01]);
                text(0,0,sprintf('stpc %1.0f',istpc));
                axis off
            end
        end
        axes('Position',[0.01,0.04,0.01,0.01]);
        text(0,0,pc_string2,'Interpreter','none');
        axis off
        axes('Position',[0.01,0.01,0.01,0.01]);
        text(0,0,tstring,'Interpreter','none');
        axis off
    end
end
