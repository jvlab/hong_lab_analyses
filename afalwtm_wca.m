function pnew=afalwtm_wca(p,w,k)
% pnew=afalwtm_wca(p,w,k) finds an equivalent affine pca model forcing 
%  the weighted sum of the affine term sum(a.*sum(w,1))=0
%
% p: p.a_uncentered, p.b_norm, p.x_centered, p.varex (frac variance explained)
%    b_norm, and x_centered are weighted PCA from data centered on columns
%    a_uncentered is the offsets required to center the data
% w: the weights
% k: the column of p.x_centered to adjust, defaults to 1
%
% pnew has all the fields of p, and also
%    pnew.a and pnew.x_true, together, which provide an equivalent model, but with
%    the weighted mean of a=0, and x_true(:,k) differeng from
%    x_centered(:,1) by a constant.
%
%   See also:  AFALTWM, AFALWTM_TEST, AFALWT_TEST.
%
if (nargin<=2) k=1; end
%
pnew=p;
%
% find how much to offset a_uncentered by
%
z=sum(p.a_uncentered.*sum(w,1))/sum(p.b_norm(k,:).*sum(w,1));
pnew.x_true=p.x_centered;
pnew.x_true(:,k)=p.x_centered(:,k)+z;
pnew.a=p.a_uncentered-z*p.b_norm(k,:);
return
