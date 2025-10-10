% afalwtm_test: test weighted affine alignment
% (multi-factor PCA with missing data)
% 
%  afalwtm allows multiple components, and an affine term.
%
%  Model:  size(d)=size(w)=[nf nr].
%   minimize sum(w.*(d-repmat(a_uncentered,nf,1)-x_centered*b_norm))^2)
%   size(a_uncentered)=[1 nr]; size(x_centered)=[nf npcs_model]; size(b_norm)=[npcs_model nr]
%
%    In afalwtm and afalwt, x_centered and b_norm are determined by
%        PCA after column-centering,
%    rows of b_norm are normalized (but need not be orthogonal for weighted pca)
%    a_uncentered (not present in afalwt) is the affine term for x_centered and b_norm
%
%   equivalent to minimizing sum(w.*(d-repmat(a,nf,1)-x_true*b_norm))^2)
%    a is the affine term for an equivalent model in which
%    sum(a.*sum(w,1))=0, and for this model, replace x_centered by x_true.
%    x_centered and x_true differ only on the first column, and only  by a
%    constant.
% 
%   the a's, b's, and x's are returned as subfields of pcastd, pcawt, etc.
%
%    See also:   AFALWT_TEST, AFALWTM_INIT, AFALWTM_ITER, AFALWTM, AFALWTM_WCA, GRMSCMDT.
%
random_seed=getinp('random seed (-1 to randomize)','d',[-1 10000],0);
if (random_seed>=0)
   rand('state',random_seed); %choose a particular seed for Matlab 5 generator
else
   rand('state',sum(100*clock)); %randomize as suggested by Matlab docs
end
%
if ~exist('nr') nr=10; end
nr=getinp('number of raters','d',[1 Inf],nr);
if ~exist('nf') nf=9000; end
nf=max(nf,nr);
nf=getinp('number of objects to rate','d',[nr Inf],nf);
if ~exist('wthr') wthr=0.8; end
wthr=getinp('0 for continuous weights, >0 for non-missing data fraction)','f',[0 1],wthr);
if ~exist('sigma') sigma=0.05; end
sigma=getinp('noise level','f',[0 Inf],sigma);
if ~exist('a_spread') a_spread=1; end
a_spread=getinp('a_spread','f',[0 Inf],a_spread);
if ~exist('b_spread') b_spread=1; end 
b_spread=getinp('b_spread','f',[0 Inf],b_spread);
if ~exist('npcs_model') npcs_model=3; end
npcs_model=getinp('npcs in model','d',[1 nr],npcs_model);
if ~exist('npcs_fit') npcs_fit=3; end
npcs_fit=getinp('npcs in fit','d',[1 10],npcs_fit);
if (npcs_fit>npcs_model)
    warning(sprintf(' number of pcs in fit (%3.0f) exceeds number of pcs in model (%3.0f)',npcs_fit,npcs_model));
end
if ~exist('eigrat') eigrat=0.75; end
eigrat=getinp('eigenvalue ratio','f',[0 1],eigrat);
if ~exist('niters') niters=5; end
niters=getinp('iterations ''by hand''','d',[0 100],niters);
if ~exist('itermax') itermax=1000; end
itermax=getinp('maximum iterations for automated routine','d',[1 Inf],itermax);
if ~exist('tol') tol=0.00001; end
tol=getinp('tolerance for automated routine','f',[0 Inf],tol);
if ~exist('ifdebug') ifdebug=0; end
ifdebug=getinp('ifdebug (for afalwtm_iter)','d',[0 1],ifdebug);
%
% make a model 
%
exact.a_uncentered=a_spread*rand(1,nr);
b=1+b_spread*rand(npcs_model,nr);
[gs,gsn]=grmscmdt(b');
b_norm=gsn';
exact.b_norm=b_norm;
clear gs gsn b b_norm;
%
% make some fake data
% d(f,r)=a_uncentered(r)+x_centered(f,m)b_norm(m,r)+gau(sigma)(f,r)
%
x=rand(nf,npcs_model);
[gs,gsn]=grmscmdt(x-repmat(mean(x,1),nf,1));
exact.x_centered=gsn.*repmat(eigrat.^[0:npcs_model-1],nf,1);
clear gs gsn
d=repmat(exact.a_uncentered,nf,1)+exact.x_centered*exact.b_norm+sigma*randn(nf,nr);
%
% do a standard pca, no data omitted
%
pcastd=afalwtm_init(d,npcs_fit,exact.b_norm);
%
% make the weight matrix
%
w=rand(nf,nr);
if (wthr>0)
    w=double(w<wthr);
    % ensure that there are at least npcs_fit ratings in each row
    for f=[find(sum(w,2)<npcs_fit)]';
        w(f,:)=double(randperm(nr)<=npcs_fit);
    end
end
%
% make the data with missing values
dm=d;
dm(find(w(:)==0))=NaN; %just to make sure we don't use it
dm_filled=dm;
%
for r=1:nr
    wnz=find(w(:,r)>0);
    dm_dot(1,r)=sum(w(wnz,r).*dm(wnz,r))/sum(w(wnz,r));
    dm_filled(find(w(:,r)<=0),r)=dm_dot(1,r);
end
clear wnz
dm_centered=dm-repmat(dm_dot,nf,1);
dm_centered_filled=dm_centered;
dm_centered_filled(find(w(:)==0))=0;
%
% do the pca by filling in missing values
pcastd_filled=afalwtm_init(dm_filled,npcs_fit,exact.b_norm);
pcastd_centered_filled=afalwtm_init(dm_centered_filled,npcs_fit,exact.b_norm);
%
% iterate a few times "by hand"
%
for k=1:niters
    if (k==1)
        [pcawt_centered{k},b_change(k)]=afalwtm_iter(pcastd_centered_filled,dm_centered,...
            w,npcs_fit,setfields([],{'iflog','ifdebug'},{1,ifdebug}));
    else
        [pcawt_centered{k},b_change(k)]=afalwtm_iter(pcawt_centered{k-1},dm_centered,...
            w,npcs_fit,setfield([],'iflog',1));
    end
    disp(pcawt_centered{k}.varex);
end
%
% do it in an automated way
%
afalwt_opts=[];
afalwt_opts=setfield(afalwt_opts,'itermax',itermax);
afalwt_opts=setfield(afalwt_opts,'tol',tol);
disp(' running afalwtm')
[pcawt,b_change,optsused]=afalwtm(dm,w,npcs_fit,afalwt_opts);
disp(optsused);
disp(' running afalwtm with noconst=1')
[pcawt_noconst,b_change_noconst,optsused_noconst]=afalwtm(dm,w,npcs_fit,setfield(afalwt_opts,'noconst',1));
disp(optsused_noconst);
