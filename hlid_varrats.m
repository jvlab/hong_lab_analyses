function vr=hlid_varrats(v)
% vr=hlid_varrats(v) computes within- and across-group variances, and their ratios
%
% v(npts,nrepts,nstims): a data array, nrepts>1, nstims>1
%
% vr.within(1,nstims): variance within repeats of each stmulus
% vr.across(1,nstims): variance between each within-repeat mean and the centroid
% vr.centroid(npts,1): the centroid, i.e., the global mean
% vr.ratio=mean(vr.across)/mean(vr.within)
% vr.frat=vr.ratio*nrepts/(1-1/nstims), distributed like F with
% vr.fdof=[(nstims-1)*npts,(nrepts-1)*nstims*npts] degreess of freedom
%
% for large nrepts,nstims, the expected value of vr.ratio is (1-1/nstims)/nrepts
%
%rng('default');c=2;r=3;s=4;N=1000000;ratio=zeros(1,N);for k=1:N;z=randn(c,r,s);ratio(k)=getfield(hlid_varrats(z),'ratio');end;F=ratio*r/(1-1/s);
%quantile(F,[0.01 .1 .5 .9 0.99]) =               0.1331    0.3493    0.9293    2.1747    4.2115
% finv([0.01 .1 0.5 .9 0.99],(s-1)*c,(r-1)*s*c) = 0.1330    0.3493    0.9300    2.1783    4.2016
% rng('default');c=5;r=2;s=3;N=1000000;ratio=zeros(1,N);for k=1:N;z=randn(c,r,s);ratio(k)=getfield(hlid_varrats(z),'ratio');end;F=ratio*r/(1-1/s);
% quantile(F,[0.01 .1 .5 .9 0.99]) =              0.2200    0.4464    0.9772    2.0595    3.8059
% finv([0.01 .1 0.5 .9 0.99],(s-1)*c,(r-1)*s*c) - 0.2194    0.4457    0.9773    2.0593    3.8049
%
% 20Mar26: add F ratio
%
%   See also:  HLID_VI_EXPLORE.
%
npts=size(v,1);
nrepts=size(v,2);
nstims=size(v,3);
%
stim_means=reshape(mean(v,2),[npts 1 nstims]);
stim_vars=sum(sum((v-repmat(stim_means,[1 nrepts 1])).^2,1),2)/(nrepts-1);
vr.within=stim_vars(:)';
%
vr.centroid=reshape(mean(stim_means,3),[npts 1]);
vr.across=sum((reshape(stim_means,[npts nstims])-repmat(vr.centroid,[1,nstims])).^2,1);
%
vr.ratio=mean(vr.across)/mean(vr.within);
vr.frat=vr.ratio*nrepts/(1-1/nstims);
vr.fdof=[(nstims-1)*npts,(nrepts-1)*nstims*npts];
return
end
