function [confmtx,aux]=xval_confmtx_make(outsamp,insamp,method)
%[confmtx,aux]=xval_confmtx_make(outsamp,insamp,method) computes a cross-validated confusion matrix 
%  with a variety of options for distances, and how to combine distances
%
% outsamp: data to be decoded: cell array of size (nstims_in,1), outsamp{istim}(isamp,idim) are the coordinates
% insamp: reference data: cell array of size (nstims_in_out,1), insamp{istim}(isamp,idim) are the coordinates
%   any row with a NaN is ignored
% method:
%   method.dist_type: distance type, 'Euclidean','1-cosine','Mahalanobis'
%   method.comb_type: how to combine across points
%      'centroid': use centroid
%      'min','max','mean','median','geomean','rms','harmonic','gravitational'
%      'powerlaw': use a power-law mean, with method.param giving the power
%      If dist_type=Mahalanobis, then comb_type must be centroid
%   method.params: parameters
%      If comb_type='powerlaw', the power law -- special cases: 1: mean, -1: harmonic mean, 2: rms, -2: gravitational
%      comb_type='powerlaw' with method.param=-Inf  is the same as comb_type='min'
%      comb_type='powerlaw' with method.param=+Inf  is the same as comb_type='max'
%      comb_type='powerlaw' with method.param=0  is the same as comb_type='geomean'
% If called with no inputs, then confmtx returns the available options
%  
% confmtx: size(nstims_out,nstims_in), confmtx(out,in) is the number of times that
%  a response to a stimulus of type out is decoded as in. sum(confmtx(out,:)) should equal size(outsamp{istim},1)
% aux.desc_full: full text description of method
% aux.desc_brief: brief text description of method
% aux.msgs: warning messages
% aux.dist_type_use,comb_type_use,params_use: methods used internally
% aux.disparity(itrial,istim_in): disparity between an out-of-sample trial
%   and an in-sample, after combining the distances according to comb_type
%   (lowest value of aux.disparity(itrial,:) yields the decode)
% aux.fc: overall fraction correct
%
%  See also:  COOTODSQ, MAHAL, XVAL_CONFIGS_MAKE.
%
dist_type_list={'Euclidean','1-cosine','Mahalanobis'};
comb_type_list={'centroid','min','max','mean','median','geomean','rms','harmonic','gravitational','powerlaw'};
params_needed=struct; %fields are number of params needed for each comb_type
params_needed.powerlaw=1;
params_default.powerlaw=2; %rms
%
special_pwr=struct;
special_pwr.mean=1;
special_pwr.rms=2;
special_pwr.harmonic=-1;
special_pwr.gravitational=-2;
%
for ic=1:length(comb_type_list)
    if ~isfield(params_needed,comb_type_list{ic})
        params_needed.(comb_type_list{ic})=0;
    end
end
aux=struct;
if (nargin<1)
    confmtx=struct;
    confmtx.dist_type_list=dist_type_list;
    confmtx.comb_type=comb_type_list;
    confmtx.params_needed=params_needed;
else
    nstims_out=size(outsamp,1);
    ndims_out=size(outsamp{1},2);
    nstims_in=size(insamp,1);
    ndims_in=size(insamp{1},2);
    %
    confmtx=zeros(nstims_out,nstims_in);
    aux.msgs=[];
    if isempty(strmatch(method.dist_type,dist_type_list,'exact'))
        aux.msgs=strvcat(aux.msgs,sprintf('distance type %s not recognized',method.dist_type));
    end
    if isempty(strmatch(method.comb_type,comb_type_list,'exact'))
        aux.msgs=strvcat(aux.msgs,sprintf('combine type %s not recognized',method.comb_type));
    end
    if ndims_out~=ndims_in
        aux.msgs=strvcat(aux.msgs,'unequal dimensions of out-of-sample and in-sample data');
    end
    if strcmp(method.dist_type,'Mahalanobis') & ~strcmp(method.comb_type,'centroid')
        aux.msgs=strvcat(aux.msgs,'Mahalanobis can only be used with centroid combination');
    end
    if ~isempty(aux.msgs)
        return
    end
    nparams=params_needed.(method.comb_type);
    if nparams>0
        method=filldefault(method,'params',params_default.(method.comb_type));
    end
    aux.desc_full=sprintf('%s, combine by %s',method.dist_type,method.comb_type);
    if nparams>0
        aux.desc_full=cat(2,aux.desc_full,sprintf(' %7.3f',method.params));
    end
    aux.desc_brief=strrep(aux.desc_full,'combine by ','');
    aux.desc_brief=strrep(aux.desc_brief,'Euclidean','Eucl');
    aux.desc_brief=strrep(aux.desc_brief,'Mahalanobis','Mahal');
    %translate to the method that will be used
    dist_type_use=method.dist_type;   
    comb_type_use=method.comb_type;
    params_use=method.params
    if strcmp(method.comb_type,'powerlaw')
        if method.params(1)==-Inf
            comb_type_use='min';
        elseif method.params(1)==Inf
            comb_type_use='max';
        elseif method.params(1)==0
            comb_type_use='geomean';
        end
    end
    if isfield(special_pwr,method.comb_type)
        comb_type_use='powerlaw';
        params_use=special_pwr.(method.comb_type);
    end
    params_use=params_use(1:params_needed.(comb_type_use));
    aux.dist_type_use=dist_type_use;
    aux.comb_type_use=comb_type_use;
    aux.params_use=params_use;
    %
    %remove all nans from insamp, ensure sufficient number of samples if Mahalanobis
    %
    insamps_count=zeros(nstims_in,1);
    for istim=1:nstims_in
        nonans=find(all(~isnan(insamp{istim}),2));
        insamp{istim}=insamp{istim}(nonans,:);
        insamps_count(istim)=length(nonans);       
    end
    if strcmp(dist_type_use,'Mahalanobis') & min(insamps_count)<=ndims_in
        aux.msgs=strvcat(aux.msgs,sprintf('Mahalanobis distance will be degenerate, ndims=%2.0f but min(nstims)=%2.0f',ndims_in,min(insamps_count)));
        return;
    end
    %
    %create concatenated outsamp as outsamp_all and tags, outsamp_tag
    %
    [outsamp_count,outsamp_all,outsamp_tag]=xval_concat_util(outsamp);
    %
    %if combine by centroid, compute centroids of insample and then disparities to insample
    %
    disparity=zeros(sum(outsamp_count),nstims_in);
    if strcmp(comb_type_use,'centroid');
        insamp_centroids=zeros(nstims_in,ndims_in);
        for istim=1:nstims_in
            insamp_centroids(istim,:)=mean(insamp{istim},1);
            if strcmp(dist_type_use,'Mahalanobis')
                disparity(:,istim)=sqrt(mahal(outsamp_all,insamp{istim}));
            end
        end
        switch dist_type_use
            case 'Euclidean'
                disparity=sqrt(cootodsq(outsamp_all,insamp_centroids));
            case '1-cosine'
                dotprods=outsamp_all*insamp_centroids';
                outsamps_normsq=sum(outsamp_all.^2,2);
                insamps_normsq=sum(insamp_centroids.^2,2);
                disparity=1-dotprods./sqrt(outsamps_normsq*insamps_normsq');
        end
    else
        %combine individual distances:  first compute the distances to each insample
        [insamp_count,insamp_all,insamp_tag]=xval_concat_util(insamp);
        switch dist_type_use
            case 'Euclidean'
                dists=sqrt(cootodsq(outsamp_all,insamp_all));
            case '1-cosine'
                dotprods=outsamp_all*insamp_all';
                outsamps_normsq=sum(outsamp_all.^2,2);
                insamps_normsq=sum(insamp_all.^2,2);
                dists=1-dotprods./sqrt(outsamps_normsq*insamps_normsq');
        end
        for istim_in=1:nstims_in
            trials_in=find(insamp_tag==istim_in);
            switch comb_type_use
                case 'min'
                    disparity(:,istim_in)=min(dists(:,trials_in),[],2);
                case 'max'
                    disparity(:,istim_in)=max(dists(:,trials_in),[],2);
                case 'median'
                    disparity(:,istim_in)=median(dists(:,trials_in),2);
                case 'geomean'
                    disparity(:,istim_in)=median(dists(:,trials_in),2);
                case 'powerlaw'
                    if params_use(1)>0
                        disparity(:,istim_in)=(mean(dists(:,trials_in).^params_use(1),2)).^(1/params_use(1));
                    else                       
                    %for negative power, any zero distance yields a zero
                        pos=find(all(dists(:,trials_in)>0,2));
                        nonpos=setdiff([1:length(trials_in)],pos);
                        disparity(nonpos,istim_in)=0;
                        disparity(pos,istim_in)=(mean(dists(pos,trials_in).^params_use(1),2)).^(1/params_use(1));                     
                    end
            end
        end
    end
    aux.disparity=disparity;
    %now decode based on minimum disparity, with fractional counts for ties
    min_disp=min(disparity,[],2);
    for istim_out=1:sum(outsamp_count)
        minlocs=find(disparity(istim_out,:)==min_disp(istim_out));
        confmtx(outsamp_tag(istim_out),minlocs)=confmtx(outsamp_tag(istim_out),minlocs)+1/length(minlocs);
    end
    %
    aux.fc=sum(diag(confmtx))/sum(confmtx(:));
end
return

function  [samp_count,samp_all,samp_tag]=xval_concat_util(samp)
%concatenate and make tags
nstims=size(samp,1);
samp_count=zeros(nstims,1);
for istim=1:nstims
    samp_count(istim)=size(samp{istim},1);
end
samp_all=zeros(sum(samp_count),size(samp{1},2));
samp_tag=zeros(sum(samp_count),1);
for istim=1:nstims
    offset=sum(samp_count([1:(istim-1)]));
    samp_tag(offset+[1:samp_count(istim)])=istim;
    samp_all(offset+[1:samp_count(istim)],:)=samp{istim};
end
return

