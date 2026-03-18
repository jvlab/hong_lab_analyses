function pv_new=hlid_vi_spatialfilter(pv,s)
% pv_new=hlid_vi_spatialfilter(pv,s) applies a spatial filter to volumetric imaging data
%
% pv: pixel values, size is [npxls time rept stim]
% s: structure created by hlid_vi_read, including
%    s.xyz_all, coordinates of the pixels
%    s.sfilt_hw: kernel half-width in each direction, 2 corresponds to [1 2 1]/4
% 
% pv_new: new pixel values
% 
%   See also:  HLID_VI_EXPLORE, HLID_VI_READ.
%
npts=size(pv,1);
nt=size(pv,2);
nrepts=size(pv,3);
nstims=size(pv,4);
%
hw=s.opts_read.sfilt_hw;
%
ker1D=binopdf([0:2*hw],2*hw,0.5);
ker=ker1D'*ker1D;
if (s.opts_read.if_log)
    disp('spatial filtering kernel:')
    disp(ker)
end
planes=unique(s.xyz_all(:,3));
pv_new=NaN(size(pv));
for iplane=1:length(planes)
    pxls_sel=find(s.xyz_all(:,3)==planes(iplane));
    xy_use=s.xyz_all(pxls_sel,[1:2]);
    if s.opts_read.if_log
        disp(sprintf(' filtering plane %3.0f: %7.0f pixels, %3.0f timepoints, %3.0f repts, %3.0f stims',...
            planes(iplane),size(xy_use,1),nt,nrepts,nstims))
    end
    xy_low=min(xy_use,[],1);
    xy_mod=xy_use-xy_low+1;
    data_mask=sparse(xy_mod(:,1),xy_mod(:,2),1);
    data_mask_full=full(data_mask);
    for it=1:nt
        for ir=1:nrepts
            for is=1:nstims
                f=sparse(xy_mod(:,1),xy_mod(:,2),pv(pxls_sel,it,ir,is));
%                f(data_mask==0)=NaN;
                fu=full(f);
                fu(data_mask_full==0)=NaN; %mask out the NaN's
                fc=conv2(fu,ker,'same');
                %put the data back
                klist=[1:length(pxls_sel)];
                for k=1:length(pxls_sel)
                    pv_new(pxls_sel(k),it,ir,is)=fc(xy_mod(k,1),xy_mod(k,2));
                end
            end
        end
    end
end
return
end
