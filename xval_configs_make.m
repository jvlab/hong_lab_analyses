function [configs,desc,opts_used]=xval_configs_make(shape,nmake,opts,defaults)
% [configs,desc,opts_used]=xval_configs_make(shape,nmake,opts,defaults) creates
%cross-validation schemes for single-trial decoding, where trials are
%parameterized by stimulus, trial number, and prep
%
% shape: shape of data array, i.e., a triple [nstims nrepts nsets]
% nmake: number of configurtions to make
%    nmake is ignored if trials are not randomly selected
% opts: options for how to set up the scheme. In some cases, NaN fields are requested from console,
%   using values in defaults, if present, as defaults
%  opts.dimnames: names of the 3 dimensions of shape, defaults to {'stim','rept','set'}
%  opts.if_log: 1 to log (default)
%  opts.max_partial_omit: maximum number of stimuli that can be held out from a single repeat
%  These options will be requested if not provided or provided as NaN
%  opts.if_single:
%     if_single(1) is not used
%     if_single(2) is 1 to force all held-out trials to be from a single repeat, otherwise 0
%     if_single(3) is 1 to force all held-out trials to be from a signle trial, otherwise 0
% defaults: default values of options  
%
% configs: array of size [nstims nrepts npreps nmake], each slice
%   (:,:,:,imake) has integers corresponding to the cross-validation 'fold' in which it is deleted
% descs: a string descriptor of the options used
% opts_used: options used
%   opts_used.[max|min]_dropped(nmake,1): [max|min]imum number of trials dropped, across all folds and sets, for configuration k (k=1:nmake)
%   opts_used.[max|min]_dropped_rept(nmake,nrpts): [max|min]imum number of trials dropped in one repeat, across all sets and folds, for configuration k (k=1:nmake)
%   opts_used.[max|min]_dropped_within(nmake,1): [max|min]imum number of trials dropped within a set, across all folds, for configuration k (k=1:nmake)
%   opts_used.[max|min]_dropped_within_rept(nmake,nrpts): [max|min]imum number of trials dropped in one repeat, across all sets and folds, for configuration k (k=1:nmake)
%
% possible future extensions:
%   * missing data
%   * more than two (or less than two) trial characterisics [now, just trial number and prep
%   * error and consistency checking
%   * choice of max_partial_omit based on dimensionality of decoding space
%
% 17Feb25: added opts_used.[max|min]_dropped*, to aid in calculation of min number of trials in any fold
% 22Feb25: added the dropped*within set specific to each set
% 10Mar25: fix bug in display of number of folds for deterministic case of one trial dropped
%
%  See also: XVAL_CONFMTX_MAKE.
%
ny_chars='NY';
%
if (nargin<=2)
    opts=struct;
end
if (nargin<=3)
    defaults=struct;
end
opts=filldefault(opts,'dimnames',{'stim','rept','set'});
opts=filldefault(opts,'if_log',1);
opts=filldefault(opts,'if_single',[0 NaN NaN]);
opts=filldefault(opts,'max_partial_omit',3);
opts=filldefault(opts,'omit_per_fold',NaN);
%
if_ok=0;
while (if_ok==0)
    desc=sprintf('[%s %s %s]=[%4.0f %4.0f %4.0f] ',opts.dimnames{:},shape);
    if opts.if_log
        disp(sprintf('creating cross-validation configurations for %s',desc));
    end
    max_partial_omit_per_fold=opts.max_partial_omit;
    for idim=2:3
        if isnan(opts.if_single(idim))
            if isfield(defaults,'if_single')
                def_val=defaults.if_single(idim);
            else
                def_val=NaN;
            end
            if isnan(def_val)
                def_val='';
            end
            opts.if_single(idim)=getinp(sprintf('1 to restrict omitted trials to a single %s',opts.dimnames{idim}),'d',[0 1],def_val);       
        end
        % if opts.if_single(idim)==0
        %     max_partial_omit_per_fold=max_partial_omit_per_fold*shape(idim);
        % end
        desc=cat(2,desc,sprintf('single %s:%s ',opts.dimnames{idim},ny_chars(opts.if_single(idim)+1)));
    end
    if isnan(opts.omit_per_fold)
        if isfield(defaults,'omit_per_fold')
            def_val=defaults.omit_per_fold;
        else
            def_val=NaN;
        end
        if isnan(def_val)
            def_val='';
        end
        opts.omit_per_fold=getinp('0 to omit all stimuli together in a fold, otherwise number of stims to omit, per repeat per set, in each fold','d',[0 max_partial_omit_per_fold],def_val);
    end
    if opts.omit_per_fold>0
        omit_string=sprintf(' omit %2.0f per %s per %s',opts.omit_per_fold,opts.dimnames{2},opts.dimnames{3});
    else
        omit_string=' omit all';
    end
    %
    if_det=0;
    if opts.omit_per_fold==0
        if_det=1;
        %deterministic case, nmake ignored
        configs=zeros([shape 1]);
        ifold=0;
        if opts.if_single(2)==1 & opts.if_single(3)==1 % one repeat, one set
            for id2=1:shape(2)
                for id3=1:shape(3)
                    ifold=ifold+1;
                    configs(:,id2,id3,1)=ifold;
                end
            end
            nfolds=ifold;
            if_ok=1;
        elseif opts.if_single(2)==1 & opts.if_single(3)==0 %one repeat, all sets
            for id2=1:shape(2)
                ifold=ifold+1;
                configs(:,id2,:,1)=ifold;
            end
            nfolds=ifold;
            if_ok=1;
        elseif opts.if_single(2)==0 & opts.if_single(3)==1 %all repeats, one set
            for id3=1:shape(3)
                ifold=ifold+1;
                configs(:,:,id3,1)=ifold;
            end
            nfolds=ifold;
            if_ok=1;
        else %opts.if_single(2)==0 & opts.if_single(3)==0 %all repeats, all sets
            disp('Disallowed setup, all trials in same fold.')
            opts.omit_per_fold=NaN;
        end
    else %a subset of stimuli held out from each trial number
        configs=zeros([shape nmake]);
        if opts.if_single(2)==1 & opts.if_single(3)==1 % one repeat, one set
            if opts.omit_per_fold==1
                ifold=prod(shape); %will allow nfolds to be computed correctly
                configs=reshape([1:prod(shape)],shape);
                if_det=1;
            else
                for imake=1:nmake
                    ifold=0;
                    for id2=1:shape(2)
                        for id3=1:shape(3)
                            seq=xval_config_util(shape(1),opts.omit_per_fold);
                            configs(:,id2,id3,imake)=seq+ifold;
                            ifold=ifold+max(seq(:));
                        end %id3
                    end %id2
                end %imake
            end
        elseif opts.if_single(2)==1 & opts.if_single(3)==0 %one repeat, all sets
            for imake=1:nmake
                ifold=0;
                for id2=1:shape(2)
                    for id3=1:shape(3)
                        seq=xval_config_util(shape(1),opts.omit_per_fold);
                        configs(:,id2,id3,imake)=seq+ifold;
                    end %id3
                    ifold=ifold+max(seq(:));
                end %id2
            end %imake
        elseif opts.if_single(2)==0 & opts.if_single(3)==1 %all repeats, one set
            for imake=1:nmake
                ifold=0;
                for id3=1:shape(3)
                    for id2=1:shape(2)
                        seq=xval_config_util(shape(1),opts.omit_per_fold);
                        configs(:,id2,id3,imake)=seq+ifold;
                    end %id2
                    ifold=ifold+max(seq(:));
                end %id3
            end %imake
        else %opts.if_single(2)==0 & opts.if_single(3)==0 %all repeats, all sets
            for imake=1:nmake
                ifold=0;
                for id3=1:shape(3)
                    for id2=1:shape(2)
                        seq=xval_config_util(shape(1),opts.omit_per_fold);
                        configs(:,id2,id3,imake)=seq+ifold;
                    end %id2
                end %id3
                ifold=ifold+max(seq(:));
            end %imake
        end
        nfolds=ifold;
        if_ok=1;
    end
end %if_ok
nconfigs_made=size(configs,4);
if (opts.if_log)
    if if_det==1
        disp(sprintf(' deterministic setup, making only one configuration'));
    end
    disp(sprintf('created %1.0f configurations, each with %4.0f folds for cross-validation',nconfigs_made,nfolds))
end
desc=cat(2,desc,sprintf('folds:%4.0f, %s',nfolds,omit_string));
nrepts=shape(2);
nsets=shape(3);
%
opts.max_dropped=zeros(nconfigs_made,1);
opts.min_dropped=zeros(nconfigs_made,1);
opts.max_dropped_rept=zeros(nconfigs_made,nrepts);
opts.min_dropped_rept=zeros(nconfigs_made,nrepts);
%
opts.max_dropped_withinset=zeros(nconfigs_made,1);
opts.min_dropped_withinset=Inf(nconfigs_made,1);
opts.max_dropped_withinset_rept=zeros(nconfigs_made,nrepts);
opts.min_dropped_withinset_rept=Inf(nconfigs_made,nrepts);
%
for imake=1:nconfigs_made
    c=configs(:,:,:,imake);
    h=hist(c(:),[1:max(c(:))]);
    opts.max_dropped(imake)=max(h); %maximum number of times a fold-tag occurs in configs(:,:,:,imake)
    opts.min_dropped(imake)=min(h(h>0)); %minimum number of times a fold-tag occurs in configs(:,:,:,imake)
    for iset=1:nsets
        c=configs(:,:,iset,imake);
        h=hist(c(:),[1:max(c(:))]);
        opts.max_dropped_withinset(imake)=max(opts.max_dropped_withinset(imake),max(h)); %maximum number of times a fold-tag occurs in configs(:,:,:,imake)
        opts.min_dropped_withinset(imake)=min(opts.min_dropped_withinset(imake),min(h(h>0))); %minimum number of times a fold-tag occurs in configs(:,:,:,imake)
    end
    for irept=1:nrepts
        cr=configs(:,irept,:,imake);
        hr=hist(cr(:),[1:max(cr(:))]);
        opts.max_dropped_rept(imake,irept)=max(hr);
        opts.min_dropped_rept(imake,irept)=min(hr(hr>0));
        for iset=1:nsets
            cr=configs(:,irept,iset,imake);
            hr=hist(cr(:),[1:max(cr(:))]);
            opts.max_dropped_withinset_rept(imake,irept)=max(opts.max_dropped_withinset_rept(imake,irept),max(hr));
            opts.min_dropped_withinset_rept(imake,irept)=min(opts.min_dropped_withinset_rept(imake,irept),min(hr(hr>0)));
        end
    end
end
opts_used=opts;
return
end

function seq=xval_config_util(nstims,omits)
%create a random sequence of length nstims, with omits or omits-1 of integers from 1 to ceil(nstims/omits)
temp=repmat([1:ceil(nstims/omits)],1,omits);
temp=temp(:); %create a sequence [1 1 1 2 2 2 3 3 3 ... ceil(nstims/omits)]'
temp=temp(randperm(length(temp)));
seq=temp([1:nstims]);
return
end


