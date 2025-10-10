function [p,b_change,optsused]=afalwtm(d,w,m,opts);
% [p,b_change,optsused]=afalwtm(d,w,m,opts) does affine alignment with missing data
%  for multiple pcs
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
% d: the data
% w: weights (same size as d)
% m: the number of principle components to fit
% opts: optional set of options
%    opts.iflog=1 to log results (defaults to 0)
%    opts.nowarnzdiv=1 to suppress warnings about zero-divides (defaults to 0)
%    opts.itermax: maximum number of iterations (defaults to 1000)
%    opts.tol:  tolerance for change in b_norm to terminate (defaults to 0.00001)
%    opts.noconst=1 to force the a-term to 0 (defaults to 0)
%
% p: p.a_uncentered, p.b_norm, p.x_centered, p.varex (frac variance explained)
%    b_norm, and x_centered are weighted PCA from data centered on columns
%    a_uncentered is the offsets required to center the data
%
%    p.a and p.x_true, together, provide an equivalent model, but with
%    the weighted mean of a=0, and x_true(:,1) differeng from
%    x_centered(:,1) by a constant.
%
% b_change: Euclidean change in b_norm on each step.
%           length(b_change)=number of iterations done
% optsused: options used
%     optsused.niters: number of iterations
%     optsused.lastchange: last change of b
%     optsused.termination='zero divide','converged','iteration limit exceeded'
%  
%   See also:  AFALWTM_TEST, AFALWTM_INIT, AFALWTM_ITER, AFALWTM_WCA, AFALWT.
%
if (nargin<=3) opts=[]; end
if ~isfield(opts,'iflog') opts.iflog=0; end
if ~isfield(opts,'nowarnzdiv') opts.nowarnzdiv=1; end
if ~isfield(opts,'itermax') opts.itermax=1000; end
if ~isfield(opts,'tol') opts.tol=0.00001; end
if ~isfield(opts,'noconst') opts.noconst=0; end
optsused=opts;
p=[];
b_change=[];
%
nf=size(d,1);
nr=size(d,2);
%
% fill unused values with the column mean, so vector operations will work
%
dfilled=d;
for r=1:nr
    wnz=find(w(:,r)==0);
    dfilled(wnz,r)=mean(d(find(w(:,r)>0),r));
end
%
if (min(sum(w,1)==0))
   optsused.termination='zero divide';
   if ~opts.nowarnzdiv
       disp(' warning: zero divide on initiation of afalwt');
   end
   return;
end
dm_dot=sum(w.*dfilled,1)./sum(w,1);
if (opts.noconst)
    dm_centered_filled=dfilled;
else
    dm_centered_filled=dfilled-repmat(dm_dot,nf,1);
end
piter=afalwtm_init(dm_centered_filled,m,ones(m,nr));
niters=0;
done=0;
while (done==0) & (niters<opts.itermax)
    niters=niters+1;
    [piter,b_change(niters),zdiv]=afalwtm_iter(piter,dm_centered_filled,w,m,opts);
    if (zdiv>0)
        optsused.termination='zero divide';
        if ~opts.nowarnzdiv
           disp(sprintf(' warning: zero divide on iteration %6.0f of afalwtm',niters));
        end
        return;
    end
    if (b_change(niters)<opts.tol)
        p=piter;
        optsused.termination='converged';
        done=1;
    end
end
if (done==0)
    p=piter;
    optsused.termination='iteration limit exceeded';
end
%
if (opts.noconst)
    p.a_uncentered=zeros(1,nr);
    p.x_true=p.x_centered;
    p.a=p.a_uncentered;
else
    p.a_uncentered=dm_dot;
    % adjust w(:,1) to center a
    p=afalwtm_wca(p,w,1);
end
optsused.niters=length(b_change);
if (length(b_change)>0)
    optsused.lastchange=b_change(end);
end
%
return
