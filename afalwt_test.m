% afalwt_test: test weighted affine alignment
% (one-factor PCA with missing data)
% 
%    See also:   AFALWT_INIT, AFALWT_ITER.
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
if ~exist('wthr') wthr=0.5; end
wthr=getinp('0 for continuous weights, >0 for non-missing data fraction)','f',[0 1],wthr);
if ~exist('sigma') sigma=0.3; end
sigma=getinp('noise level','f',[0 Inf],sigma);
if ~exist('a_spread') a_spread=1; end
a_spread=getinp('a_spread','f',[0 Inf],a_spread);
if ~exist('b_spread') b_spread=1; end
b_spread=getinp('b_spread','f',[0 Inf],b_spread);
if ~exist('niters') niters=5; end
niters=getinp('iterations ''by hand''','d',[0 100],niters);
if ~exist('itermax') itermax=1000; end
itermax=getinp('maximum iterations for automated routine','d',[1 Inf],itermax);
if ~exist('tol') tol=0.00001; end
tol=getinp('tolerance for automated routine','f',[0 Inf],tol);
%
% make a model 
%
a=a_spread*rand(1,nr);
a=a-mean(a); %center a
exact.a=a;
b=1+b_spread*rand(1,nr);
b_norm=b/sqrt(sum(b.^2));
exact.b_norm=b_norm;
%
% make some fake data
% d(f,r)=a(r)+x_true(f)*b_norm(r)+gau(sigma)(f,r)
%
exact.x_true=rand(nf,1);
d=repmat(a,nf,1)+exact.x_true*b_norm+sigma*randn(nf,nr);
%
% do a standard pca, no data omitted
%
pcastd=afalwt_init(d,exact.b_norm);
%
% make the weight matrix
%
w=rand(nf,nr);
if (wthr>0)
    w=double(w<wthr);
    % ensure that there is at least one rating in each row
    for f=[find(sum(w,2)==0)]';
        w(f,:)=double(randperm(nr)==1);
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
pcastd_filled=afalwt_init(dm_filled,exact.b_norm);
pcastd_centered_filled=afalwt_init(dm_centered_filled,exact.b_norm);
%
% iterate a few times "by hand"
%
for k=1:niters
    if (k==1)
        [pcawt_centered{k},b_change(k)]=afalwt_iter(pcastd_centered_filled,dm_centered,w,setfield([],'iflog',1));
    else
        [pcawt_centered{k},b_change(k)]=afalwt_iter(pcawt_centered{k-1},dm_centered,w,setfield([],'iflog',1));
    end
end
%
% do it in an automated way
%
afalwt_opts=[];
afalwt_opts=setfield(afalwt_opts,'itermax',itermax);
afalwt_opts=setfield(afalwt_opts,'tol',tol);
disp(' running afalwt')
[pcawt,b_change,optsused]=afalwt(dm,w,afalwt_opts);
disp(optsused);
disp(' running afalwt with noconst=1')
[pcawt_noconst,b_change_noconst,optsused_noconst]=afalwt(dm,w,setfield(afalwt_opts,'noconst',1));
disp(optsused_noconst);
