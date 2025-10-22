function pcastd=afalwt_init(d,b_norm_dir)
% pcastd=afalwt_init(d,b_norm_dir) initiates an affine alignment with missing data.
% This does a standard PCA to find the vectors a, b_norm, and x such that
% d(f,r)=a_r+b_norm_r*(x_true(f))+gau(sigma)(f,r) 
% restricted to have sum(a_r)=0 and sum(b_r^2)=1
%
%  if b_norm_dir is given, then pcastd.b_norm is multiplied by +/-1 to align with it.
% answers and auxiliary values returned as fields of pcastd
%
%   See also:  AFALWT_TEST, AFALWT_ITER, AFALWT.
%

if (nargin<=1) b_norm_dir=0; end
nf=size(d,1);
nr=size(d,2);
d_centered=d-repmat(mean(d,1),nf,1);
[u,s,v]=svd(d_centered,0); %typically nf>nr - the zero argument reduces the svd size.
s=diag(s(1:nr,1:nr));
pcastd.varex=(s(1)^2)/sum(s.^2); %variance explained by first principal component
pcastd.b_norm=v(:,1)'/sqrt(sum(v(:,1).^2));

%
%absorb the multiplier into x

%
pcastd.x_centered=u(:,1)*s(1); % First principal component. It's mean is zero because the input data had mean zero.
%disp(sprintf('boo %d, bee %d, been %d',mean(pcastd.x_centered),mean(u(:,1)),s(1))) %%%% Want to remove this eventually.
pcastd.x_true=pcastd.x_centered+repmat(mean(d(:)),nf,1)/mean(pcastd.b_norm);
pcastd.a=mean(d,1)-mean(pcastd.x_true)*pcastd.b_norm;
%
% align the sign, if a sign vector is given
%
if sum(pcastd.b_norm.*b_norm_dir)<0
    pcastd.b_norm=-pcastd.b_norm;
    pcastd.x_centered=-pcastd.x_centered;
end
%
return
