function [resps,coords_all]=hlid_fill_merge_svd(s,opts)
% [resps,f]=hlid_fill_merge_svd(s,opts):  use missing-data routine to fill in and merge ORN data
%
%  This is essentially the modularized version of processing in hild_orn_merge, followed by optional mean subtraction, and hlid_coords_svd.
%
%  s: data structure.  s{ifile}.(opts.response_name).(opts.response_name2) is the array of responses (stimuli x glomeruli)
%  opts: options
%
%  resps: array of size(length(opts.nstims_use),length(opts.nglomeruli)), responses filled in
%  coords_all: svd of coords and auxiliary data (see hlid_coords_svd)
%
%   See also:  HLID_PRED_MAGNIF_DEMO, AFALWT, HLID_ORN_MERGE.
%
nfiles_use=length(opts.files_use);
nstims=length(opts.stimuli_use);
nglomeruli=length(opts.glomeruli);
%
resps_all=zeros(nstims,nglomeruli,nfiles_use);
%
glomeruli_missing=zeros(nglomeruli,nfiles_use);
for ifile_ptr=1:nfiles_use
    ifile=opts.files_use(ifile_ptr);
    glomeruli_missing(opts.nancols{ifile},ifile_ptr)=1;
    resps_all(:,:,ifile_ptr)=s{ifile}.(opts.response_name).(opts.response_name2)(opts.stimuli_use,:);
end
for ipresent=0:nfiles_use
    if opts.if_log
        disp(sprintf('number of glomeruli present in at least %2.0f datasets: %3.0f (missing in %2.0f or less)',...
            ipresent,sum(sum(glomeruli_missing,2)<=(nfiles_use-ipresent)),nfiles_use-ipresent));
    end
end
glomeruli_use=find(sum(glomeruli_missing,2)<=nfiles_use-opts.min_present);
nglomeruli_use=length(glomeruli_use);
resps_gu=resps_all(:,glomeruli_use,:);
%
%fill in missing data by affine interpolation
% 
resps_gur=reshape(resps_gu,[nstims*nglomeruli_use,nfiles_use]);
if ~exist('afalwt_opts') afalwt_opts=struct;end
resps_tofill=isnan(resps_gur);
if_canfill=1;
if any(all(resps_tofill==1,1))
    if opts.if_log
        disp('cannot fill in missing data, no stimuli present for some (glomerulus,file) pair')
    end
    if_canfill=0;
end
if any(all(resps_tofill==1,2))
    if opts.if_log
        disp('cannot fill in missing data, no (glomerulus,file) pair present for some stimulus')
    end
    if_canfill=0;
end
if if_canfill
    [afalwt_fit,afalwt_b_change,afalwt_optsused]=afalwt(resps_gur,1-resps_tofill,afalwt_opts);
    resps_gur_fitted=(afalwt_fit.x_true*afalwt_fit.b_norm+repmat(afalwt_fit.a,size(resps_gur,1),1)); %interpolated data
    resps_gur_filled=resps_gur;
    resps_gur_filled(resps_tofill)=resps_gur_fitted(resps_tofill);
    resps_gu_filled=reshape(resps_gur_filled,[nstims nglomeruli_use nfiles_use]);
else
    resps_gu_filled=resps_gu;
end
%
resps=reshape(afalwt_fit.x_true,[nstims nglomeruli_use]); %use regression slope as response measure
if opts.if_restore_size
    resps=resps*geomean(afalwt_fit.b_norm);;
end
if (opts.if_submean)
    resps=resps-repmat(mean(resps,1),nstims,1);
end
%
if opts.if_log
    disp(sprintf('%4.0f NaN values filled in',sum(resps_tofill(:))));
end
%
%do SVD
%
maxdim_allowed=min(size(resps))-opts.if_submean;
maxdim_use=maxdim_allowed;
[fnew,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(struct(),resps,maxdim_allowed,maxdim_use,opts.if_submean,[],[],setfield(struct(),'if_log',opts.if_log));
if opts.if_plot
    %
    %show data, before glomerulus selection and after selection, normalized within each dataset
    %
    for ifig=1:4
        switch ifig
            case 1
                figname='all raw data';
                resps_plot=resps_all;
                glomeruli_plot=[1:nglomeruli];
                nsubp=nfiles_use;
            case 2
                figname='raw data from selected glomeruli';
                resps_plot=resps_gu;
                glomeruli_plot=glomeruli_use;
                nsubp=nfiles_use;
            case 3
                figname='raw data from selected glomeruli with missing data filled in';
                resps_plot=resps_gu_filled;
                glomeruli_plot=glomeruli_use;
                nsubp=nfiles_use;
            case 4
                figname='merged response measure';
                resps_plot=resps;
                glomeruli_plot=glomeruli_use;
                nsubp=1;
        end
        figure;
        set(gcf,'Position',[50 100 1800 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figname);
        [nr,nc]=nicesubp(nsubp);
        for ifile_ptr=1:nsubp
            ifile=opts.files_use(ifile_ptr);
            subplot(nr,nc,ifile_ptr);
            imagesc(resps_plot(:,:,ifile_ptr));
            minmax=[min(min(resps_plot(:,:,ifile_ptr),[],'omitnan')),max(max(resps_plot(:,:,ifile_ptr),[],'omitnan'))];
            title_string=sprintf('%s  [%6.3f %6.3f]',strrep(opts.filenames_short{ifile},'.mat',''),minmax);
            title(title_string,'Interpreter','none');
            set(gca,'FontSize',7);
            set(gca,'XTick',[1:length(glomeruli_plot)]);
            set(gca,'XTickLabel',opts.glomeruli(glomeruli_plot));
            set(gca,'YTick',[1:nstims]);
            set(gca,'YTickLabel',opts.stim_labels(opts.stimuli_use));
        end
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,figname,'Interpreter','none');
        axis off;
    end %ifig
end


return
end
