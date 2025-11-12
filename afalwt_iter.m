function [pnew,b_change,zdiv]=afalwt_iter(p,d,w,opts)
% [pnew,b_change,zdiv]=afalwt_iter(p,d,w,opts) does the iteration step of w weighted affine alignment
% with missing data (no constant term) -- i.e., missing-data pca with first PC only
%
% d(f,r)=b_norm_r*(x_true(f))+gau(sigma)(f,r), weighted by w(f,r) 
%  if w(f,r)=0, then d(f,r) can be NaN
%
%  uses projection pursuit method
%    Shum, H-Y., Ikeuchi, K., and Reddy, R. (1995)
%   Principal component analysis with missing data and its appliation to
%    polyhedral object modeling.  IEEE Trans. Pattern Analysis and Machine
%    Intelligence 17, 854-867.
%
% p: guess from previous cycle, with fields p.b_norm and p.x_true
% d: the data.  It should be centered (i.e., mean weighted by w on each column is 0)
% w: weights (same size as d)
% opts: optional set of options
%    opts.iflog=1 to log results (defaults to 0)
%    opts.nowarnzdiv=1 to suppress warnings about zero-divides (defaults to 0)
%
% pnew: next guess, with fields pnew.b_norm, pnew.x_true, pnew.varex (frac variance explained)
% b_change: Euclidean change in b_norm
% zdiv:  1 if there was a zero-divide exception (not enough data to do fit)
%
% 12Nov25: bug fix:  add opts.legacy, defaults to 1, set to 0 to actually use Euclidean norm (per Jon Drover)
%
%   See also:  AFALWT_TEST, AFALWT_INIT, AFALWT.
%
if (nargin<=3) opts=[]; end
if ~isfield(opts,'iflog') opts.iflog=0; end
if ~isfield(opts,'nowarnzdiv') opts.nowarnzdiv=0; end
if ~isfield(opts,'legacy') opts.legacy=1; end
%
nf=size(d,1);
nr=size(d,2);
%
dfilled=d;
dfilled(find(w(:)==0))=0;
zdiv=0;
%
% calculate new x from old b
%
x_den=sum(w.*repmat(p.b_norm.^2,nf,1),2);
if (min(x_den)==0)
    if (opts.nowarnzdiv==0)
        disp(' warning: x_den=0 at indices')
        disp(find(x_den==0)');
    end
    zdiv=1;
end
x=sum(w.*dfilled.*repmat(p.b_norm,nf,1),2)./x_den;
%
% calculate new b from new x
%
b_den=sum(w.*repmat(x.^2,1,nr),1);
if (min(b_den)==0)
    if (opts.nowarnzdiv==0)
        disp(' warning: b_den=0 at indices')
        disp(find(b_den==0));
    end
    zdiv=1;
end
b=sum(w.*dfilled.*repmat(x,1,nr),1)./b_den;
%
%absorb the multiplier into x
%
pnew.b_norm=b./sqrt(sum(b.^2));
pnew.x_true=x.*sqrt(sum(b.^2));
%
% align the sign to previous
%
if sum(pnew.b_norm.*p.b_norm)<0
    pnew.b_norm=-pnew.b_norm;
    pnew.x_true=-pnew.x_true;
end
d_pred=pnew.x_true*pnew.b_norm;
pnew.varex=1-sum(sum(w.*(d_pred-dfilled).^2))/sum(sum(w.*(dfilled.^2)));
%
if opts.legacy %added 12Nov25
    b_change=sqrt(sum(pnew.b_norm-p.b_norm).^2);
else
    b_change=sqrt(sum((pnew.b_norm-p.b_norm).^2));
end
if (opts.iflog)
    disp(sprintf(' Euclidean change in pnew.b_norm is %10.8f:',b_change))
    disp(pnew.b_norm);
end
return
