function [pnew,b_change,zdiv]=afalwtm_iter(p,d,w,m,opts)
% [pnew,b_change,zdiv]=afalwtm_iter(p,d,w,m,opts) does the iteration step of w weighted affine alignment
% with missing data (no constant term) -- i.e., missing-data pca, many pc's
%
% d(f,r)=a_r+{sum over m} b_norm_(r)<m>*(x_true(f)<m>)+gau(sigma)(f,r) 
%  if w(f,r)=0, then d(f,r) can be NaN
%
%  uses projection pursuit method
%    Shum, H-Y., Ikeuchi, K., and Reddy, R. (1995)
%   Principal component analysis with missing data and its appliation to
%    polyhedral object modeling.  IEEE Trans. Pattern Analysis and Machine
%    Intelligence 17, 854-867.
%
% p: guess from previous cycle, with fields p.b_norm and p.x_centered (not like afalwt_iter)
%    size(p.b_norm)    =[m size(d,2)];
%    size(p.x_centered)=[size(d,1) m];
% d: the data.  It should be centered (i.e., mean weighted by w on each column is 0)
% w: weights (same size as d)
% m: the number of principle components to fit
% opts: optional set of options
%    opts.iflog=1 to log results (defaults to 0)
%    opts.nowarnzdiv=1 to suppress warnings about zero-divides (defaults to 0)
%
% pnew: next guess, with fields pnew.b_norm, pnew.x_centered, pnew.varex (frac variance explained)
% b_change: Euclidean change in b_norm
% zdiv:  1 if there was a zero-divide exception (not enough data to do fit)
%
%   See also:  AFALWTM_TEST, AFALWTM_INIT, AFALWT_ITER.
%
if (nargin<=4) opts=[]; end
if ~isfield(opts,'iflog') opts.iflog=0; end
if ~isfield(opts,'ifdebug') opts.ifdebug=0; end
if ~isfield(opts,'nowarnzdiv') opts.nowarnzdiv=0; end
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
for g=1:nf
    wb=repmat(w(g,:),m,1).*p.b_norm;
    Bprime=p.b_norm*wb'; %B=wb*p.b_norm';
    wd=w(g,:).*dfilled(g,:);
    Aprime=p.b_norm*wd'; %A=wd*p.b_norm'
    x(g,:)=(Bprime\Aprime)'; %x(g,:)=(B'\A')';
end
if (opts.ifdebug)
    if (m==1)
        %do 1-d case, like afalwt_iter.m
        x_den=sum(w.*repmat(p.b_norm.^2,nf,1),2);
        x_1d=sum(w.*dfilled.*repmat(p.b_norm,nf,1),2)./x_den;
        disp('   prev      new     new(1d)');
        disp([p.x_centered x x_1d]);
    else
        disp('   prev    new');
        disp([p.x_centered x]);
    end
end
%
% calculate new b from new x
%
for s=1:nr
    wx=repmat(w(:,s),1,m).*x;
    Dprime=x'*wx; % D=wx'*x;
    wd=w(:,s).*dfilled(:,s);
    Cprime= x'*wd; % C=wd'*x;
    b(:,s)=(Dprime\Cprime);% b(:,s)=D'\C';
end
if (opts.ifdebug)
    if (m==1)
        b_den=sum(w.*repmat(x.^2,1,nr),1);
        b_1d=sum(w.*dfilled.*repmat(x,1,nr),1)./b_den;
        disp('  prev''     new''    new(1d)''');
        disp([p.b_norm' b' b_1d']);
    else
        disp('  prev''   new''');
        disp([p.b_norm' b']);
    end
end
%
%absorb the multiplier into x
%
renorm=sqrt(sum(b.^2,2));
pnew.b_norm=b./repmat(renorm,1,nr);
pnew.x_centered=x.*repmat(renorm',nf,1);
%
% align the sign to previous
%
for k=1:m
    if sum(pnew.b_norm(k,:).*p.b_norm(k,:))<0
        pnew.b_norm(k,:)=-pnew.b_norm(k,:);
        pnew.x_centered(:,k)=-pnew.x_centered(:,k);
    end
end
for k=1:m %calculate cumulative variance explained, and then amount for each
    %since no guarantee of orthogonality if weights are nontrivial
    d_pred=pnew.x_centered(:,1:k)*pnew.b_norm(1:k,:);
    var_unex=sum(sum(w.*(d_pred-dfilled).^2))/sum(sum(w.*(dfilled.^2)));
    if (k==1)
        pnew.varex(k)=1-var_unex;
    else
        pnew.varex(k)=1-sum(pnew.varex(1:k-1))-var_unex;
    end
end
%
b_change=sqrt(sum(pnew.b_norm(:)-p.b_norm(:)).^2);
if (opts.iflog)
    disp(sprintf(' Euclidean change in pnew.b_norm is %10.8f:',b_change))
    disp(pnew.b_norm);
end
return
