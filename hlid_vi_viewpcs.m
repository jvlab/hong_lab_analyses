%hlid_vi_viewpcs: view principal components of KC volumetric imaging, data from George Barnum, Hong Lab
%
%   See also:  HLID_VI_READ, HLID_VI_SPATIALFILTER, HLID_VARRATS, HLID_VI_PREPROC.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='gbarnum_mb247_soma_20241027_a_test_1.hdf5'; end
if ~exist('nstims') nstims=24; end
if ~exist('nrepts') nrepts=5; end
if ~exist('n_tbins') n_tbins=5; end %default number of time bins to use
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
rept_list=getinp('repeat list','d',[1 nrepts],1:nrepts);
stim_list=getinp('stimulus list','d',[1 nstims],1:nstims);
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
        bin_ranges=[(1+[0 cumsum(tbins(1:end-1))]);cumsum(tbins)];
        spatem=reshape(svd_u(:,ipc),[n_pixels resp_minlength]);
        repstm=reshape(svd_v(:,ipc),[n_repts n_stims]);
        spatem_binned=zeros(n_pixels,n_tbins);
        for tbin=1:n_tbins
            spatem_binned(:,tbin)=mean(spatem(:,[bin_ranges(1,tbin):bin_ranges(2,tbin)]),2);
        end
        %plot with a plane in each row
        figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,pc_string,' ',tstring));
        set(gcf,'Position',[50 50 1400 800]);
        spatem_binned_max=max(abs(spatem_binned(:)));
        for iplane=1:n_planes
            %code borrowed from hlid_vi_spatialfilter
            pxls_sel=find(xyz(:,3)==planes(iplane));
            xy_use=xyz(pxls_sel,[1:2]);
            xy_min=min(xy_use,[],1);
            xy_max=max(xy_use,[],1);
            xy_mod=xy_use-xy_min+1;
            data_mask=sparse(xy_mod(:,1),xy_mod(:,2),1);
            data_mask_full=full(data_mask);
            for tbin=1:n_tbins
                f=sparse(xy_mod(:,1),xy_mod(:,2),spatem_binned(pxls_sel,tbin));
                fu=full(f);
                fu(data_mask_full==0)=NaN; %mask out the NaN's
                subplot(n_planes,n_tbins,tbin+(iplane-1)*n_tbins);
                switch if_hv
                    case 1
                        if if_unifscale
                            imagesc(fu,spatem_binned_max*[-1 1]);
                        else
                            imagesc(fu,max(abs(fu(:)))*[-1 1]);
                        end
                        xlabel(sprintf('fr [%2.0f %2.0f]',bin_ranges(:,tbin)));
                        ylabel(sprintf('plane %2.0f',planes(iplane)));
                        axis equal;
                        set(gca,'YTick',1+[0 xy_max(1)-xy_min(1)]);
                        set(gca,'YTickLabels',[xy_min(1),xy_max(1)]);
                        set(gca,'XTick',1+[0 xy_max(2)-xy_min(2)]);
                        set(gca,'XTickLabels',[xy_min(2),xy_max(2)]);
                        axis tight;
                    case 2
                        if if_unifscale
                            imagesc(fu,spatem_binned_max*[-1 1]);
                        else
                            imagesc(fu,max(abs(fu(:)))*[-1 1]);
                        end
                        xlabel(sprintf('fr [%2.0f %2.0f]',bin_ranges(:,tbin)));
                        ylabel(sprintf('plane %2.0f',planes(iplane)));
                        axis equal;
                        set(gca,'XTick',1+[0 xy_max(1)-xy_min(1)]);
                        set(gca,'XTickLabels',[xy_min(1),xy_max(1)]);
                        set(gca,'YTick',1+[0 xy_max(2)-xy_min(2)]);
                        set(gca,'YTickLabels',[xy_min(2),xy_max(2)]);
                        axis tight;
                end %if_hv
            end %tbin
        end %iplane
        axes('Position',[0.01,0.04,0.01,0.01]);
        text(0,0,pc_string2,'Interpreter','none');
        axis off
        axes('Position',[0.01,0.01,0.01,0.01]);
        text(0,0,tstring,'Interpreter','none');
        axis off
    end
end
