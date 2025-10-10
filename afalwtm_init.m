function pcastd=afalwtm_init(d,m,b_norm_dir)
% pcastd=afalwtm_init(d,m,b_norm_dir) initiates an affine alignment with missing data.
% This does a standard PCA to find the vectors a, b_norm, and x such that
%
% d(f,r)=a_r+{sum over m} b_norm_(r)<m>*(x_true(f)<m>)+gau(sigma)(f,r) 
% 
% sum(b_norm(r)^2)=1
%
% here, c_centered is calculated after weiighted centering of the data
% but in afalwt_init, a is zero-centered after weighing (sum(a_r*w(r,f))=0)
% so this does not reduce to afalwt if m=1.
%
% d assumed to be tall and narrow
% m is number of eigenvalues sought
%
% for each row of b_norm_dir that is given, pcastd.b_norm is multiplied by +/-1 to align with it.
% answers and auxiliary values returned as fields of pcastd
%
%   See also:  AFALWTM_TEST.
%
if (nargin<=2) b_norm_dir=0; end
nf=size(d,1);
nr=size(d,2);
pcastd.a_uncentered=mean(d,1);
d_centered=d-repmat(pcastd.a_uncentered,nf,1);
[u,s,v]=svd(d_centered,0); %typically nf>nr
s=diag(s(1:m,1:m));
pcastd.varex=(s(1:m).^2)/sum(s.^2); %variance explained by each principal component
pcastd.b_norm=(v(:,1:m)./repmat(sqrt(sum(v(:,1:m).^2,1)),nr,1))';
%
%absorb the multiplier into x
pcastd.x_centered=u(:,1:m).*repmat(s(1:m)',nf,1);
%
% align
for k=1:m
    if sum(pcastd.b_norm(k,:).*b_norm_dir(k,:))<0
        pcastd.b_norm(k,:)=-pcastd.b_norm(k,:);
        pcastd.x_centered(:,k)=-pcastd.x_centered(:,k);
    end
end
return
