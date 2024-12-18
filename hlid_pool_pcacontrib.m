%hlid_pool_pcacontrib: analyze the contributions of each dataset to a pooled dataset
% 
%   See also:  HLID_RASTIM2COORDS_POOL.
% 
if ~exist('pool_fn_def') pool_fn_def='./data/kc_soma_nls/hlid_odor17_coords_kc_soma_nls_pooled.mat'; end
fn=getinp('file name','s',[],pool_fn_def);
f=load(fn);
aux=f.coord_opts.aux;
nsets=length(aux.nresps_each);
npcs=size(aux.v,2);
disp(aux);
p=zeros(nsets,sum(aux.nresps_each));
for iset=1:nsets
    p(iset,sum(aux.nresps_each(1:(iset-1)))+[1:aux.nresps_each(iset)])=1/aux.nresps_each(iset);
    disp(sprintf(' set %1.0f: %30s, %4.0f responses in pool',iset,f.dsid(iset,:),aux.nresps_each(iset)));
end
wt_mean=p*aux.v;
wt_rms=sqrt(p*(aux.v.^2));
figure;
set(gcf,'Position',[100 100 800 600]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','contribs to pca');
%mean
subplot(2,1,1);
plot(wt_mean');
legend(f.dsid,'Location','Best','FontSize',7);
set(gca,'XLim',[0 npcs]);
set(gca,'XTick',[1:npcs]);
xlabel('pc')
title('mean weight');
%rms
subplot(2,1,2);
plot(wt_rms');
legend(f.dsid,'Location','Best','FontSize',7);
set(gca,'XLim',[0 npcs]);
set(gca,'XTick',[1:npcs]);
xlabel('pc')
title('rms weight');


