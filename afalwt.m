function [p,b_change,optsused]=afalwt(d,w,opts)
% [p,b_change,optsused]=afalwt(d,w,opts) does affine alignment with missing data
%
% d(f,r)=a_r+b_norm_r*(x_true(f))+gau(sigma)(f,r), weighted by w(f,r) 
%  if w(f,r)=0, then d(f,r) can be NaN
%
%  uses projection pursuit method
%    Shum, H-Y., Ikeuchi, K., and Reddy, R. (1995)
%   Principal component analysis with missing data and its application to
%    polyhedral object modeling.  IEEE Trans. Pattern Analysis and Machine
%    Intelligence 17, 854-867.
%
% d: the data
% w: weights (same size as d)
% opts: optional set of options
%    opts.iflog=1 to log results (defaults to 0)
%    opts.nowarnzdiv=1 to suppress warnings about zero-divides (defaults to 0)
%    opts.itermax: maximum number of iterations (defaults to 1000)
%    opts.tol:  tolerance for change in b_norm to terminate (defaults to 0.00001)
%    opts.noconst=1 to force the a-term to 0 (defaults to 0)
%
% p: p.b_norm, p.x_true, p.varex (frac variance explained)
% b_change: Euclidean change in b_norm on each step.
%           length(b_change)=number of iterations done
% optsused: options used
%     optsused.niters: number of iterations
%     optsused.lastchange: last change of b
%     optsused.termination='zero divide','converged','iteration limit exceeded'
%  
%   See also:  AFALWT_TEST, AFALWT_INIT, AFALWT_ITER, MTC_SOID_XLSRUN_SUMM, FIGGND_DBASE_SUMM.
%
if (nargin<=2) opts=[]; end
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
   if ~opts.nowarnzdiv % This seems like a reason to stop.
       disp(' warning: zero divide on initiation of afalwt');
   end
   return;
end
%%% Weighted mean, by default we want to remove this.
dm_dot=sum(w.*dfilled,1)./sum(w,1);
if (opts.noconst)
    dm_centered_filled=dfilled;
else
    dm_centered_filled=dfilled-repmat(dm_dot,nf,1);
end
%%%%%% INIT CALL %%%%%%%
piter=afalwt_init(dm_centered_filled,ones(1,nr));
niters=0;
done=0;
while (done==0) & (niters<opts.itermax)
    niters=niters+1;
    %%%%%%% ITER CALL %%%%%%%%
    [piter,b_change(niters),zdiv]=afalwt_iter(piter,dm_centered_filled,w,opts);
    if (zdiv>0)
        optsused.termination='zero divide';
        if ~opts.nowarnzdiv
           disp(sprintf(' warning: zero divide on iteration %6.0f of afalwt',niters));
        end
        return; %% This returns, but doesn't terminate. Should it? Does it terminate in the caller?
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
p.x_centered=p.x_true; %calculation was done on centered values
if (opts.noconst)
    p.a=zeros(1,nr);
else
    %x_dot=(sum(sum(w.*dfilled))./sum(sum(w)))/mean(p.b_norm);
    x_dot=sum(sum(w.*dfilled))./sum(sum(w.*repmat(p.b_norm,nf,1)));
    p.x_true=p.x_centered+x_dot;
    p.a=dm_dot-p.b_norm*x_dot;
end
optsused.niters=length(b_change);
if (length(b_change)>0)
    optsused.lastchange=b_change(end);
end
%
return
