%hlid_rastim_trial_decode: single-trial decoding analysis
%
% does a cross-validated decode, beginning with response space construction with trials left out
% uses consensus space to align data across preps
% 
% Constructs a representational space using single-trial responses,
% optionally subtracting mean, optionally normalizing the magnitude,
% (setting Euclidean length to 1, so that distances are monotonically related to 1-correlation)
%
% if_singleprep:
%   0 (default): decodes by aligning across preps
%   1: decodes by using the repeats within the same prep, and does not attempt to align across preps
% if_noembed:
%   0 (default): only does decoding based on embedding
%   1: also does decoding based on ROI data, without a dimension-reduction, requires if_singleprep=1
% if_embedbyprep:
%   0 (default):  embedding does an SVD on each repeat separately, and then aligns them (nrepts_gp=1)
%   1: embedding is carried out all repeats of the same prep together. (nrepts_gp=nrepts).
%    With this option, alignment only will align responses to the same stimulus on the same repeat.
% If both if_singleprep=1 and if_embedbyprep=1, then only some kinds of cross-validation configurations are allowed,
%    to prevent all responses to the same stimulus being dropped for the in-sample of a fold.
% 
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO,
%  PROCRUSTES_CONSENSUS, PROCRUSTES_COMPAT, HLID_RASTIM_TRIAL_VIS,
%  XVAL_CONFIGS_MAKE, XVAL_CONFMTX_MAKE, HLID_RASTIM_TRIAL_READ,
%  HLID_RASTIM_TRIAL_DECODE_3DPLOT, HLID_RASTIM_TRIAL_DECODE_CONFMTX.
%
% 15Feb25: add a single-prep mode (decode only within each prep): if_singleprep
% 20Feb25: allow for dimension list to be non-contiguous
% 03Mar25: begin to add option for single-prep decoding without embedding, if_noembed: dimlist entry=Inf
% 05Mar25: now defaults to 'econ' svd (can change by setting if_econ_svd=0)
% 09Mar25: begin option to do embedding of all repeats of a prep at the same time: if_embedbyprep
% 10Mar25: fix bug when all stimuli in a trial are dropped
% 14Mar25: finished option to do embedding of all repeats of a prep at the same time: if_embedbyprep
% 17Mar25: better logic to allow if_singleprep=1 with if_embedbyprep=1
%
hlid_setup;  %invoke hlid_localopts; set up opts_read and opts_plot
if_debug=getinp('1 for debug mode','d',[0 1]);
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if_singleprep=getinp('1 for decoding only within single preps','d',[0 1],0);
if if_singleprep==1
    if_noembed=getinp('1 to also decode without embedding','d',[0 1],0);
else
    if_noembed=0;
end
if_embedbyprep=getinp('1 to embed all repeats of a prep together','d',[0 1],0);
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
if ~exist('sub_labels') 
    if if_debug
        sub_labels={''};
    else
        sub_labels={'',' (mean sub)'}; end %subtract mean from responses, can replace by a subset to shorten analysis
end
nsubs=length(sub_labels);
if ~exist('preproc_labels') 
    if if_debug
        preproc_labels={'raw'};
    else
        preproc_labels={'raw','normalized'}; end %can replace by a subset to shorten analysis
end
npreprocs=length(preproc_labels);
if ~exist('if_econ_svd')
    if_econ_svd=1;
end
%
if ~exist('dec_methods')
    dec_methods=cell(0);
    dec_methods{1}.dist_type='Euclidean';
    dec_methods{1}.comb_type='mean';
    dec_methods{2}.dist_type='Euclidean';
    dec_methods{2}.comb_type='centroid';
    dec_methods{3}.dist_type='1-cosine';
    dec_methods{3}.comb_type='centroid';
    if if_singleprep %single-fly: Mahalanobis will fail because there are not enough trials
        dec_methods{4}.dist_type='Euclidean';
        dec_methods{4}.comb_type='gravitational';
    else
        dec_methods{4}.dist_type='Mahalanobis';
        dec_methods{4}.comb_type='centroid';
    end
end
ndecs=length(dec_methods);
dec_labels=cell(1,ndecs);
%
if ~exist('opts_pcon') opts_pcon=struct; end %consensus options
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
if ~exist('opts_xv') 
    opts_xv=struct; 
    opts_xv=filldefault(opts_xv,'if_single',[0 NaN NaN]);
    if if_singleprep
        opts_xv.if_single(3)=1; %restrict to within single prep
    end
    if if_embedbyprep
        opts_xv.if_single(2)=1; %restrict to a single repeat (irrelevant, only one repeat, so don't even ask)
    end
end %for cross-validation
if ~exist('xv_defaults') xv_defaults=struct; end
if ~exist('xv_nmake') xv_nmake=1; end
%
if ~exist('nrepts') nrepts=3; end %number of repeats
if ~exist('nsets_min') nsets_min=1; end
if ~exist('dmax') 
    if if_debug
        dmax=3;
    else
        dmax=7;%max representational space to create
    end
end
if ~exist('dmin') 
    dmin=1;
end
%
if ~exist('if_log') if_log=0; end
hlid_rastim_trial_read;
%
% metadata=cell(nfiles,1);
% dsids=cell(nfiles,1);
% resps_mean=cell(nfiles,1); mean responses within stimuli(nstims,nrois_avail)
% trial_ptrs=cell(nfiles,1); array of (nstims,nrepts)
% resps_trial=cell(nfiles,1); trial-by-trial responses, stimuli unscrambled, (nstims,nrepts,nrois_avail)
% trial_sequence=cell(nfiles,1); stimulus sequence, stimuli as strings
% stims_avail=cell(nfiles,1); list of available stimuli in each file, beginning at 1
% rois_avail=cell(nfiles,1); list of roi numbers kept for analysis, beginning at 1
% rois=cell(nfiles,1):  original rois
% nrois_avail(ifile): number of rois available
%
for k=1:nstims
    stimulus_names_display{k}=stimulus_names(k,1:-1+min(find(stimulus_names(k,:)==' ')));
end
%
disp(sprintf('number of files: %4.0f',nfiles));
if if_singleprep
    nsets=1;
else
    nsets=getinp('number of datasets to analyze together','d',[nsets_min nfiles]);
end
%
subsamps_list=nchoosek(1:nfiles,nsets);
nsubsamps=size(subsamps_list,1);
if (if_singleprep)
    nsubsamps_use=getinp('number of subsample sets to use','d',[1 nsubsamps],nsubsamps);
    subsamps_use=[1:nsubsamps];
else
    nsubsamps_use=getinp('number of subsample sets to use','d',[1 nsubsamps],1);
    subsamps_use=randperm(nsubsamps);
end
if_all_subsamps=double(nsubsamps_use==nsubsamps);
subsamps_use=subsamps_use(1:nsubsamps_use);
%
if if_embedbyprep==0
    nrepts_gp=1;
else
    nrepts_gp=nrepts;
end
max_embed=nstims*nrepts_gp;
%
%logic to get size of reprentational space (embedding dimension), and
%cross-validation schemes, and check for consistency -- max embed dimension
%cannot exceed number of responses retained in any insample
%
if_ok=0;
while (if_ok==0)
    dmax=getinp('max dimension of representational space to create','d',[max(1,dmin) max_embed],dmax);
    dmin=getinp('min dimension of representational space to create','d',[1 dmax],dmin);
    dimlist=getinp('list of dimensions to create','d',[dmin dmax],[dmin:dmax]);
    xv_nmake=getinp('number of cross-validation configurations to make','d',[1 Inf],xv_nmake);
    ixv_valid=0;
    while ixv_valid==0
        if if_embedbyprep==0
            [xv_configs,xv_label,opts_xv_used]=xval_configs_make([nstims,nrepts,nsets],xv_nmake,opts_xv,xv_defaults);
            xv_configs_embed=xv_configs; %embedding is by rept, not by prep
        else
            opts_xv_embedbyprep=opts_xv;
            opts_xv_embedbyprep.blocks_allowed=[1 nstims]; %either delete one stimulus or an entire repeat
            opts_xv_embedbyprep.phases_allowed=0; %one starting phase only
            [xv_configs_embed,xv_label_temp,opts_xv_used]=xval_configs_make([nstims*nrepts,1,nsets],xv_nmake,opts_xv_embedbyprep,xv_defaults);
            xv_configs=reshape(xv_configs_embed,[nstims nrepts nsets size(xv_configs_embed,4)]);
            xv_label=strrep(xv_label_temp,'[stim rept set]','[stimXrept 1 set]');
        end
        %verify that in each cross-validation configuration, no stimulus is completely dropped
        if any(any(min(min(xv_configs,[],2),[],3)==max(max(xv_configs,[],2),[],3)))
            disp('configuration is not valid, some stimulus is completely dropped in at least one cross-validation set')
        else
            ixv_valid=1;
        end
    end
    xv_nmake=size(xv_configs,4); %in case this was a deterministic configuration, and only one xv_config was made
    %
    %determine if any embeddings of in-sample data are missing all examples of (stim x rept)
    %
    xv_configs_exclude=cell(xv_nmake,max(xv_configs(:)));
    xv_configs_exclude_count=zeros(xv_nmake,max(xv_configs(:)));
    for ixv_make=1:xv_nmake
        which_drop=xv_configs_embed(:,:,:,ixv_make);
        folds_alldropped=unique(which_drop(find(min(min(which_drop,[],2),[],3)==max(max(which_drop,[],2),[],3)))); %folds that have all examples of (stim x rept) dropped
        for fold_ptr=1:length(folds_alldropped)
            ifold=folds_alldropped(fold_ptr);
            xv_configs_exclude{ixv_make,ifold}=find(which_drop==ifold); %  xv_configs_exclude{ixv_make,ifold} is a list of the stim x rept that are completely missing in any fold
                           
            xv_configs_exclude_count(ixv_make,ifold)=sum(which_drop==ifold);
        end
        if sum(xv_configs_exclude_count(ixv_make,:))>0
            disp(sprintf('Note that  for xv config %2.0f, %2.0f folds drop all examples of stimuli x repeats',ixv_make,sum(xv_configs_exclude_count>0)));
            disp(sprintf('     fold: %s',sprintf(' %3.0f',folds_alldropped)));
            disp(sprintf(' ndropped: %s',sprintf(' %3.0f',xv_configs_exclude_count(ixv_make,folds_alldropped))));
        end
    end
    max_dropped=max(opts_xv_used.max_dropped_withinset_rept(:));
    min_dropped=min(opts_xv_used.min_dropped_withinset_rept(:));
    if (max_dropped==max_embed) & (min_dropped==max_embed) 
        max_dropped=0;
    end %if all stimuli of a repeat are dropped together, then no fold has fewer trials
    dmax_allowed=max_embed-max_dropped;
    disp(sprintf('Some insample pca analyses will have as few as %2.0f trials (%2.0f dropped).  Max rep space dimension is %1.0f.',...
        dmax_allowed,max_dropped,dmax));
    if dmax_allowed>=dmax
        if_ok=1;
    else
        disp('Need to adjust cross-validation or max rep space dimension.')
        dmin=1;
    end
end
if if_noembed==1
    dimlist=[dimlist Inf];
end
%
opts_pcon.allow_scale=getinp('1 to allow scaling','d',[0 1],1);
scaling_token=(opts_pcon.allow_scale==1);
if opts_pcon.allow_reflection==1
    reflection_token='best';
else
    reflection_token=false;
end
%
if_3dplot=getinp('1 to plot response space for 3d fit, last subsamp, last fold, last xv set','d',[0 1]);
%
%do the preprocessing on all files, as this is independent of the
%sub-selection each subsample
%
resps_alltrials=cell(nsubs,npreprocs,nfiles);
for ifile=1:nfiles
    for isub=1:nsubs
        rs=resps_mean{ifile};
        rt=resps_trial{ifile};
        switch sub_labels{isub}
            case ''
            case ' (mean sub)'
            rs_xm=mean(rs,1,'omitnan'); %global mean
            rt=rt-repmat(reshape(rs_xm,[1 1 nrois_avail(ifile)]),[nstims nrepts 1]);
        end
        for ipreproc=1:npreprocs
            switch preproc_labels{ipreproc}
                case 'raw'
                case 'normalized'
                    rt_norm=sqrt(sum(rt.^2,3));
                    rt_norm(rt_norm==0)=1;
                    rt=rt./repmat(rt_norm,[1 1 nrois_avail(ifile)]);
            end
            resps_alltrials{isub,ipreproc,ifile}=rt;
        end %ipreproc
    end %isub
    disp(sprintf('preprocessed file %s',dsids{ifile}));
end %ifile
clear resps_mean resps_trial rois rs rs_xm rt rt_norm
%initialize confusion matrices
confusion_matrices=zeros(nstims,nstims,ndecs,length(dimlist),nsubs,npreprocs,nsubsamps_use); % d1: actual stim, d2: decoded stim, d3: decision rule, d4: id_ptr, d5: sub mean d6: normalize, d7: subsample set
maxerr_insamp_recon=repmat(-Inf,[nsubs npreprocs nsubsamps_use]); %maximum reconstruction error when number of pcs=number of in-sample responses
for isubsamp=1:nsubsamps_use
    subsamp_sel=subsamps_list(subsamps_use(isubsamp),:);
    disp(sprintf(' subsample %3.0f of %3.0f: original datasets %s',isubsamp,nsubsamps_use,sprintf(' %3.0f ',subsamp_sel)));
    nrois=nrois_avail(subsamp_sel(:)); %number of rois for each set selected
    for isub=1:nsubs %mean subtract?
        for ipreproc=1:npreprocs
            for id_ptr=1:length(dimlist)
                id=dimlist(id_ptr);
                if_pca=and(id>=dmin,id<=dmax);
                if if_pca
                    npcs=id;
                    %
                    %create a consensus from all of the data, to be used just for initialization
                    %
                    coords_nodrop=zeros(nstims*nrepts_gp,id,nrepts*nsets/nrepts_gp); %new for embedbyprep
                    for iset=1:nsets
                        resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                        for irept=1:nrepts/nrepts_gp
                            irepts=irept+[0:nrepts_gp-1];
                            rpca=reshape(resps(:,irepts,:),nstims*nrepts_gp,nrois_avail(subsamp_sel(iset)));
                            nonans_pca=find(all(~isnan(rpca),2));
                            if if_econ_svd
                                [u_nonan,s,v]=svd(rpca(nonans_pca,:),'econ'); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                            else
                                [u_nonan,s,v]=svd(rpca(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                            end
                            u=nan(size(rpca,1),npcs);
                            u(nonans_pca,:)=u_nonan(:,1:npcs);
                            % coords_nodrop(:,:,irept+(iset-1)*nrepts)=u*s(1:npcs,1:npcs);
                            coords_nodrop(:,:,irept+(iset-1)*nrepts/nrepts_gp)=u*s(1:npcs,1:npcs);
                        end
                    end
                    consensus_nodrop=procrustes_consensus(coords_nodrop,opts_pcon);
                    disp(sprintf('no-drop consensus [%2.0f x %2.0f] made for isubsamp %2.0f  %12s %10s dim %2.0f',...
                        size(consensus_nodrop,1),size(consensus_nodrop,2),isubsamp,sub_labels{isub},preproc_labels{ipreproc},id));
                else
                    disp(sprintf('no-embed decoding      for isubsamp %2.0f  %12s %10s',isubsamp,sub_labels{isub},preproc_labels{ipreproc}));
                end %if_pca
                %now do each fold:
                % eliminate some stimuli from conesnsus calculation
                % use consensus_nodrop to initialize
                % retain the transforms so that the dropped stimuli can be decoded
                %
                for ixv_make=1:xv_nmake
                    nfolds=max(max(max(xv_configs(:,:,:,ixv_make))));
                    for ifold=1:nfolds
                        % if_pca=1: for embedding with pca, do in-sample first and then out-of-sample, since one needs to transform the out-of-sample
                        % data by the consensus alignment based on in-sample 
                        % if_pca=0: no embedding, single-prep only:  do in-sample and out-of-sample together
                        ndrop_thisfold=sum(sum(sum(double(xv_configs(:,:,:,ixv_make)==ifold)))); %number dropped on this fold
                        coords_insample_align=cell(nstims,1); %coords_insample_align{istim}(itrial,id) are in-sample coords of istim
                        if if_pca
                            %coords_insample=zeros(nstims,id,nrepts*nsets-ndrop_thisfold); % eliminated 03Mar25
                            coords_insample=zeros(nstims*nrepts_gp,id,nrepts*nsets/nrepts_gp); %out-of-sample data will be replaced by NaN's after PCA
                            trans_lookup=zeros(nrepts/nrepts_gp,nsets);
                            itrial_ct=0;
                            v_insample=cell(nrepts/nrepts_gp,nsets);
                            droplist=cell(nrepts/nrepts_gp,nsets); %what is dropped
                            %omit the dropped stimuli and do pca to
                            %create an in-sample space
                            trials_havedata=[];
                            itrial_ct=0;
                            for iset=1:nsets
                                resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                                 for irept=1:nrepts/nrepts_gp
                                    irepts=irept+[0:nrepts_gp-1];
                                    itrial_ct=itrial_ct+1;
                                    rpca=reshape(resps(:,irepts,:),nstims*nrepts_gp,nrois(iset));
                                    droplist{irept,iset}=find(xv_configs(:,irepts,iset,ixv_make)==ifold); %which stimuli are dropped from this rept and set; OK if > irepts is a vector
                                    rpca(droplist{irept,iset},:)=NaN;
                                    nonans_pca=find(all(~isnan(rpca),2));
                                     if any(~isnan(rpca(:)))
                                        if if_econ_svd
                                            [u_nonan,s,v]=svd(rpca(nonans_pca,:),'econ'); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                                        else
                                            [u_nonan,s,v]=svd(rpca(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                                        end
                                        u=nan(size(rpca,1),npcs);
                                        u(nonans_pca,:)=u_nonan(:,1:npcs);
                                        %
                                        v_insample{irept,iset}=v(:,[1:npcs]); %for out-of-sample, coords=resp*v
                                        trials_havedata=[trials_havedata itrial_ct];
                                        coords_insample(:,:,itrial_ct)=u*s(1:npcs,1:npcs);                           
                                        trans_lookup(irept,iset)=itrial_ct; %where to find the transform
                                        %
                                        if npcs==length(nonans_pca)
                                            recon=coords_insample(:,:,itrial_ct)*v(:,1:npcs)';
                                            recon_error=max(max(abs(recon(nonans_pca,:)-rpca(nonans_pca,:))));
                                            maxerr_insamp_recon(isub,ipreproc,isubsamp)=max(maxerr_insamp_recon(isub,ipreproc,isubsamp),recon_error);
                                        end
                                    end
                                 end %irept
                            end %iset
                            ntrials_tot=itrial_ct;
                            %find consensus of every insample set of trials, and the transforms to each ts_insample{trans_lookup(irept,iset)}
                            %after excluding the stims (or stims x repts) that are not present
                            %  xv_configs_exclude{ixv_make,ifold} is a list of the stim x rept that are completely missing in any fold
                            stims_keep=setdiff([1:size(coords_insample,1)],xv_configs_exclude{ixv_make,ifold});
                            opts_pcon.initial_guess=consensus_nodrop(stims_keep,:);
                            consensus_insample=NaN(nstims*nrepts_gp,npcs);
                            znew_insample=NaN(nstims*nrepts_gp,npcs,length(trials_havedata));
                            [consensus_insample(stims_keep,:),znew_insample(stims_keep,:,:),ts_insample,details_insample,opts_pcon_insample_used]=...
                                procrustes_consensus(coords_insample(stims_keep,:,trials_havedata),setfield(opts_pcon,'initialize_set',0));
                            %
                            for istim=1:nstims
                                znew_stim=[];
                                for irept=1:nrepts_gp %stack up the responses to the same stimulus on different uns
                                    znew_stim=[znew_stim;reshape(znew_insample(istim+(irept-1)*nstims,:,:),[id,length(trials_havedata)])']; %d1: itrial, d2: dim
                                end
%                                znew_stim=reshape(znew_insample(istim,:,:),id,ntrials_tot)'; %d1: itrial, d2: dim
                                which_trials=find(all(~isnan(znew_stim),2)); %a subset of [1:nsets*nrepts]
                                coords_insample_align{istim}=znew_stim(which_trials,:); %d1: trial (excluding those dropped), d2: dim
                            end %end
                            %
                            %transform the dropped stimuli into the in-sample activity space and decode
                            %
                            % if only some of a trial are dropped, use the transformation from the in-sample 
                            % if the entire trial is dropped, align to the in-sample consensus
                            %trans_lookup(irept,iset) will be zero if the entire trial is omitted; otherwise, it points to
                            %ts_insample{*} that transforms that trial into consensus_insample
                            %
                            coords_outsample_align=cell(nstims,1); %coords_outsample_align{istim}(itrial,id) are out-of-sample coords of istim
                            for iset=1:nsets
                                resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                                for irept=1:nrepts/nrepts_gp;
                                    irepts=irept+[0:nrepts_gp-1];
                                    rpca=reshape(resps(:,irepts,:),nstims*nrepts_gp,nrois(iset));
                                    ndrop=length(droplist{irept,iset});
                                    %if only some of the trial is dropped, then use the alignment from the rest of the trial
                                    if ndrop>0 & ndrop<=opts_xv_used.omit_per_fold
                                        resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                                        %rpca_outsample=reshape(resps(droplist{irept,iset},irept,:),[ndrop,nrois_avail(subsamp_sel(iset))]);
                                        resps_gp=reshape(resps,[nstims*nrepts_gp,nrepts/nrepts_gp,nrois(iset)]);
                                        rpca_outsample=reshape(resps_gp(droplist{irept,iset},irept,:),[ndrop,nrois(iset)]);
                                        coords_outsample=rpca_outsample*v_insample{irept,iset}; %unaligned coords
                                        coords_outsample_align_alldrop=psg_geomodels_apply('procrustes',coords_outsample,...
                                            procrustes_compat(ts_insample{trans_lookup(irept,iset)})); %d1: each dropped stim in droplist, d2: coord
                                        %put these into coords_outsample_alignment
                                        for idrop=1:ndrop
                                            istim=droplist{irept,iset}(idrop);
                                            coords_outsample_align{mod(istim-1,nstims)+1}(end+1,:)=coords_outsample_align_alldrop(idrop,:);
                                        end
                                        %if_embedbyprep mods to here
                                    elseif ndrop>opts_xv_used.omit_per_fold %if all of the trial is dropped (corresponds to omit_per_fold=0)
                                        %transform entire repeat into consensus_insample as a reference, first doing a private pca
                                        rpca_outsample=reshape(resps(:,irepts,:),[nstims*nrepts_gp,nrois(iset)]);
                                        nonans_pca=find(all(~isnan(rpca_outsample),2));
                                        if if_econ_svd
                                            [u_nonan,s,v]=svd(rpca_outsample(nonans_pca,:),'econ'); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                                        else
                                            [u_nonan,s,v]=svd(rpca_outsample(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                                        end
                                        u=nan(size(rpca,1),npcs);
                                        u(nonans_pca,:)=u_nonan(:,1:npcs);
                                        coords_outsample=u*s(1:npcs,1:npcs);
                                        %align with consensus_insample, using
                                        %an option for procrustes that matches the option for procrustes_consensus
                                        coords_outsample_xform=NaN(nstims,npcs);
                                        [d_procrust,coords_outsample_xform(nonans_pca,:),xform]=...
                                            procrustes(consensus_insample(nonans_pca,:),coords_outsample(nonans_pca,:),'Scaling',scaling_token,'Reflection',reflection_token);
                                        if (opts_pcon.allow_offset==0) %remove offset if requested
                                            coords_outsample_xform(nonans_pca,:)=coords_outsample_xform(nonans_pca,:)-xform.c;
                                        end
                                        for idrop=1:ndrop
                                            istim=droplist{irept,iset}(idrop);
                                            coords_outsample_align{mod(istim-1,nstims)+1}(end+1,:)=coords_outsample_xform(istim,:);
                                        end
                                    end % ndrop
                                end %irept
                            end %iset
                        else %if_pca=0, no embedding. nsets=1
                            iset=1;
                            resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)};
                            coords_insample=NaN(nstims,nrois(iset),nrepts*nsets);
                            coords_outsample=NaN(nstims,nrois(iset),nrepts*nsets);
                            coords_insample_align=cell(nstims,1); %coords_insample_align{istim}(itrial,id) are in-sample coords of istim
                            coords_outsample_align=cell(nstims,1); %coords_outsample_align{istim}(itrial,id) are out-of-sample coords of istim
                            for irept=1:nrepts %get the responses for the dropped nd non-dropped stimuli on each repeat
                                droplist{irept,iset}=find(xv_configs(:,irept,iset,ixv_make)==ifold); %which stimuli are dropped from this rept and set
                                keeplist=setdiff([1:nstims],droplist{irept,iset}); %which stimuli are kept
                                coords_insample(keeplist,:,irept)=reshape(resps(keeplist,irept,:),[length(keeplist),size(resps,3)]); %d1: stim, d2: roi, d3: rept
                                coords_outsample(droplist{irept,iset},:,irept)=reshape(resps(droplist{irept,iset},irept,:),[length(droplist{irept,iset}),size(resps,3)]); %d1: stim, d2: roi, d3: rept
                            end
                            for istim=1:nstims
                                znew_stim=reshape(coords_insample(istim,:,:),[size(resps,3) nrepts])'; % d1: itrial, d2: roi
                                which_trials=find(all(~isnan(znew_stim),2));
                                coords_insample_align{istim}=znew_stim(which_trials,:);
                                %
                                znew_stim_out=reshape(coords_outsample(istim,:,:),[size(resps,3) nrepts])';
                                which_trials_out=find(all(~isnan(znew_stim_out),2));
                                coords_outsample_align{istim}=znew_stim_out(which_trials_out,:);
                            end
                        end %if_pca
                        %
                        % add to confusion matrices                
                        % confusion_matrices=zeros(nstims,nstims,ndecs,length(timlist),nsubs,npreprocs,nsubsamps_use); % d1: actual stim, d2: decoded stim, d3: decision rule, d4: id_ptr, d5: sub mean d6: normalize, d7: subsample set
                        for idec=1:ndecs
                            [confmtx,aux]=xval_confmtx_make(coords_outsample_align,coords_insample_align,dec_methods{idec});
                            dec_labels{idec}=aux.desc_brief;
                            confusion_matrices(:,:,idec,id_ptr,ipreproc,isub,isubsamp)=confusion_matrices(:,:,idec,id_ptr,ipreproc,isub,isubsamp)+confmtx;
                        end
                    end %ifold
                end %ixv_make
                if (if_3dplot) & (id==3)
                    hlid_rastim_trial_decode_3dplot;
                end
            end %id_ptr
        end %ipreproc
    end %isub (mean subtract)   
    if any(maxerr_insamp_recon(:)~=-Inf)
        disp('maximum reconstruction error when number of pcs is number of in-sample trials')
        disp(maxerr_insamp_recon(:,:,isubsamp));
    end
end %isubsamp
%
%save results
%
results=struct;
results.nrepts=nrepts;
results.nstims=nstims;
results.nfiles=nfiles;
results.nsets=nsets;
results.nsubs=nsubs;
results.sub_labels=sub_labels;
results.npreprocs=npreprocs;
results.preproc_labels=preproc_labels;
results.dec_methods=dec_methods;
results.dec_labels=dec_labels;
results.ndecs=ndecs;
results.dmax=dmax;
results.dmin=dmin;
results.dimlist=dimlist;
results.if_frozen=if_frozen;
results.if_debug=if_debug;
results.if_singleprep=if_singleprep;
results.if_noembed=if_noembed;
results.if_embedbyprep=if_embedbyprep;
results.if_all_subsamps=if_all_subsamps;
results.if_econ_svd=if_econ_svd;
%
results.metadata=metadata;
results.dsids=dsids;
results.stims_avail=stims_avail; %list of available stimuli in each file, beginning at 1
results.rois_avail=rois_avail;
results.stimulus_names=stimulus_names;
results.stimulus_names_display=stimulus_names_display;
%
results.subsamps_list_used=subsamps_list(subsamps_use,:);
%
results.xv_configs=xv_configs;
results.xv_label=xv_label;
results.opts_xv=opts_xv;
results.opts_xv_used=opts_xv_used;
results.opts_pcon=opts_pcon;
%
results.confusion_matrices=confusion_matrices;
results.confusion_matrices_dims={'d1: actual stim, d2: decoded stim, d3: decision rule, d4: id_ptr, d5: sub mean d6: normalize, d7: subsample set'};
%
disp('results structure created, consider saving it');
%
disp('confusion matrices can be plotted with hlid_rastim_trial_decode_confmtx')
