function [p,b_change,optsused]=afalwtmg(d,w,m,opts);
% [p,b_change,optsused]=afalwtmg(d,w,m,opts) does affine alignment with missing data
%  for multiple pcs via a greedy method
%
%  It finds one PC at a time, subtracting each off, and then re-analyzing
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
% b_change: Euclidean change in b_norm for last extracted component
%           length(b_change)=number of iterations done
% optsused: options used
%     optsused.niters: number of iterations
%     optsused.niters_greedy: number of iterations at each step
%     optsused.lastchange: last change of b
%     optsused.termination='zero divide','converged','iteration limit exceeded'
%     optsused.termination_greedy{:}:  how terminated on each extracted component
%     optsused.b_change_greedy{:}: Euclidean change in b_norm for each extracted component
%  
%   See also:  AFALWTM, AFALWTMG_TEST, AFALWTM_INIT, AFALWTM_ITER, AFALWTM_WCA, AFALWT.
%
if (nargin<=3) opts=[]; end
if ~isfield(opts,'iflog') opts.iflog=0; end
if ~isfield(opts,'nowarnzdiv') opts.nowarnzdiv=1; end
if ~isfield(opts,'itermax') opts.itermax=1000; end
if ~isfield(opts,'tol') opts.tol=0.00001; end
if ~isfield(opts,'noconst') opts.noconst=0; end
optsused=opts;
%
nf=size(d,1);
nr=size(d,2);
%
d_step=d;
if (opts.iflog==1)
    disp(sprintf('afalwtmg greedy algorithm step %3.0f',1));
end
[pstep,b_change,optsused]=afalwtm(d_step,w,1,optsused);
p=pstep;
optsused.termination_greedy{1}=optsused.termination;
optsused.niters_greedy=optsused.niters;
optsused.b_change_greedy{1}=b_change;
%
if (m==1)
    return;
end
% now, remove each pc found, and look for the next
im=1;
for mg=2:m
    if (strmatch(optsused.termination,'converged'))
        im=mg; %one more convergence
        if (opts.iflog==1)
            disp(sprintf('afalwtmg greedy algorithm step %3.0f',mg));
        end
        %subtract off most recent component
        %p.a_uncentered, p.b_norm, p.x_centered
        d_step=d_step-repmat(pstep.a_uncentered,nf,1)-pstep.x_centered*pstep.b_norm;
        [pstep,b_change,optsused]=afalwtm(d_step,w,1,optsused);
        optsused.termination_greedy{mg}=optsused.termination;
        optsused.b_change_greedy{mg}=b_change;
        optsused.niters_greedy(mg)=optsused.niters;
        p.b_norm(mg,:)=pstep.b_norm;
    	p.x_centered(:,mg)=pstep.x_centered;
        p.x_true(:,mg)=pstep.x_true;
        p.a=p.a+pstep.a; %these components add
        p.a_uncentered=p.a_uncentered+pstep.a_uncentered;
        p.varex(mg)=pstep.varex*(1-sum(p.varex)); %fraction of variance as yet unexplained
    else
        if (opts.iflog==1)
            disp(sprintf('afalwtmg greedy algorithm step %3.0f skipped because of prior nonconvergence',mg));
        end
    end
end
%
if (im<m)
    for mg=(im+1):m
        optsused.termination_greedy{mg}='nonconvergence at earlier stage';
        optsused.b_change_greedy{mg}=[];
        optsused.niters_greedy(mg)=0;
    end
    p.b_norm((im+1):m,:)=0;
	p.x_centered(:,(im+1):m)=0;
    p.x_true(:,(im+1):m)=0;
    p.varex((im+1):m)=0;
end
return
