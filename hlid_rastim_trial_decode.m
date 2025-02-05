%hlid_rastim_trial_decode: single-trial decoding analysis
%
% does a cross-validated decode, beginning with response space construction with trials left out
% uses consensus space to align data across individuals
% 
% Constructs a representational space using single-trial or stimulus-averaged raw (z-scored) responses,
% optionally subtracting mean, optionally normalizing the magnitude,
% (setting Euclidean length to 1, so that distances are monotonically related to 1-correlation),
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO,
%  PROCRUSTES_CONSENSUS, PROCRUSTES_COMPAT, HLID_RASTIM_TRIAL_VIS,
%  XVAL_CONFIGS_MAKE, XVAL_CONFMTX_MAKE,
%  HLID_RASTIM_TRIAL_DECODE_3DPLOT, HLID_RASTIM_TRIAL_DECODE_CONFMTX.
%
hlid_setup;  %invoke hlid_localopts; set up opts_read and opts_plot
if_debug=getinp('1 for debug mode','d',[0 1]);
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
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
        sub_labels={'',' (mean sub)'}; end %can replace by a subset to shorten analysis
end
nsubs=length(sub_labels);
if ~exist('preproc_labels') 
    if if_debug
        preproc_labels={'raw'};
    else
        preproc_labels={'raw','normalized'}; end %can replace by a subset to shorten analysis
end
npreprocs=length(preproc_labels);
%
if ~exist('dec_methods')
    dec_methods=cell(0);
    dec_methods{1}.dist_type='Euclidean';
    dec_methods{1}.comb_type='mean';
    dec_methods{2}.dist_type='Euclidean';
    dec_methods{2}.comb_type='centroid';
    dec_methods{3}.dist_type='1-cosine';
    dec_methods{3}.comb_type='centroid';
    dec_methods{4}.dist_type='Mahalanobis';
    dec_methods{4}.comb_type='centroid';
end
ndecs=length(dec_methods);
dec_labels=cell(1,ndecs);
%
if ~exist('opts_pcon') opts_pcon=struct; end %consensus options
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
if ~exist('opts_xv') opts_xv=struct; end %for cross-validation
if ~exist('xv_defaults') xv_defaults=struct; end
if ~exist('xv_nmake') xv_nmake=1; end
%
if ~exist('nrepts') nrepts=3; end %number of repeats
if ~exist('nsets_min') nsets_min=2; end
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
nsets=getinp('number of datasets to analyze together','d',[nsets_min nfiles]);
%
subsamps_list=nchoosek(1:nfiles,nsets);
nsubsamps=size(subsamps_list,1);
nsubsamps_use=getinp('number of subsample sets to use','d',[1 nsubsamps],1);
subsamps_use=randperm(nsubsamps);
subsamps_use=subsamps_use(1:nsubsamps_use);
%
dmax=getinp('max dimension of representational space to create','d',[max(1,dmin) nstims],dmax);
dmin=getinp('min dimension of representational space to create','d',[1 dmax],dmin);
xv_nmake=getinp('number of cross-validation configurations to make','d',[1 Inf],xv_nmake);
[xv_configs,xv_label,opts_xv_used]=xval_configs_make([nstims,nrepts,nsets],xv_nmake,opts_xv,xv_defaults);
%
opts_pcon.allow_scale=getinp('allow scaling','d',[0 1],1);
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
confusion_matrices=zeros(nstims,nstims,ndecs,dmax,nsubs,npreprocs,nsubsamps_use); % d1: actual stim, d2: decoded stim, d3: decision rule, d4: dmax, d5: sub mean d6: normalize, d7: subsample set
for isubsamp=1:nsubsamps_use
    subsamp_sel=subsamps_list(subsamps_use(isubsamp),:);
    disp(sprintf(' subsample %3.0f of %3.0f: original datasets %s',isubsamp,nsubsamps_use,sprintf(' %3.0f ',subsamp_sel)));
    for isub=1:nsubs %mean subtract?
        for ipreproc=1:npreprocs
            for id=dmin:dmax
                npcs=id;
                %
                %create a consensus from all of the data, to be used just for initialization
                %
                coords_nodrop=zeros(nstims,id,nrepts*nsets); 
                for iset=1:nsets
                    resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                    for irept=1:nrepts
                        rpca=reshape(resps(:,irept,:),nstims,nrois_avail(subsamp_sel(iset)));
                        nonans_pca=find(all(~isnan(rpca),2));
                        [u_nonan,s,v]=svd(rpca(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                        u=nan(size(rpca,1),npcs);
                        u(nonans_pca,:)=u_nonan(:,1:npcs);
                        coords_nodrop(:,:,irept+(iset-1)*nrepts)=u*s(1:npcs,1:npcs);
                    end
                end
                consensus_nodrop=procrustes_consensus(coords_nodrop,opts_pcon);
                disp(sprintf('no-drop consensus made for isubsamp %2.0f  %12s %10s dim %2.0f',isubsamp,sub_labels{isub},preproc_labels{ipreproc},id));
                %now do each fold:
                % eliminate some stimuli from conesnsus calculation
                % use consensus_nodrop to initialize
                % retain the transforms so that the dropped stimuli can be decoded
                for ixv_make=1:xv_nmake
                    nfolds=max(max(max(xv_configs(:,:,:,ixv_make))));
                    for ifold=1:nfolds
                        ndrop_thisfold=sum(sum(sum(double(xv_configs(:,:,:,ixv_make)==ifold)))); %number dropped on this fold
                        coords_insample=zeros(nstims,id,nrepts*nsets-ndrop_thisfold);
                        trans_lookup=zeros(nrepts,nsets);
                        itrial_ct=0;
                        v_insample=cell(nrepts,nsets);
                        droplist=cell(nrepts,nsets); %what is dropped
                        %omit the dropped stimuli and do pca to
                        %create an in-sample space
                        for iset=1:nsets
                            resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                             for irept=1:nrepts
                                rpca=reshape(resps(:,irept,:),nstims,nrois_avail(subsamp_sel(iset)));
                                droplist{irept,iset}=find(xv_configs(:,irept,iset,ixv_make)==ifold); %which stimuli are dropped from this rept and set
                                rpca(droplist{irept,iset},:)=NaN;
                                nonans_pca=find(all(~isnan(rpca),2));
                                 if any(~isnan(rpca(:)))
                                    [u_nonan,s,v]=svd(rpca(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
                                    u=nan(size(rpca,1),npcs);
                                    u(nonans_pca,:)=u_nonan(:,1:npcs);
                                    %
                                    v_insample{irept,iset}=v(:,[1:npcs]); %for out-of-sample, coords=resp*v
                                    %
                                    itrial_ct=itrial_ct+1;
                                    coords_insample(:,:,itrial_ct)=u*s(1:npcs,1:npcs);
                                    trans_lookup(irept,iset)=itrial_ct; %where to find the transform
                                end
                             end %irept
                        end %iset
                        ntrials_tot=itrial_ct;
                        opts_pcon.initial_guess=consensus_nodrop;
                        %find consensus of every insample set of trials,
                        %and the transforms to each ts_insample{trans_lookup(irept,iset)}
                        [consensus_insample,znew_insample,ts_insample,details_insample,opts_pcon_insample_used]=...
                            procrustes_consensus(coords_insample,setfield(opts_pcon,'initialize_set',0));
                        %
                        coords_insample_align=cell(nstims,1); %coords_insample_align{istim}(itrial,id) are in-sample coords of istim
                        for istim=1:nstims
                            znew_stim=reshape(znew_insample(istim,:,:),id,ntrials_tot)'; %d1: itrial, d2: dim
                            which_trials=find(all(~isnan(znew_stim),2));
                            coords_insample_align{istim}=znew_stim(which_trials,:);
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
                            for irept=1:nrepts
                                rpca=reshape(resps(:,irept,:),nstims,nrois_avail(subsamp_sel(iset)));
                                ndrop=length(droplist{irept,iset});
                                %if only some of the trial is dropped, then use the alignment from the rest of the trial
                                if ndrop>0 & ndrop<=opts_xv_used.omit_per_fold
                                    resps=resps_alltrials{isub,ipreproc,subsamp_sel(iset)}; %stim, irept, iroi
                                    rpca_outsample=reshape(resps(droplist{irept,iset},irept,:),[ndrop,nrois_avail(subsamp_sel(iset))]);
                                    coords_outsample=rpca_outsample*v_insample{irept,iset}; %unaligned coords
                                    coords_outsample_align_alldrop=psg_geomodels_apply('procrustes',coords_outsample,...
                                        procrustes_compat(ts_insample{trans_lookup(irept,iset)})); %d1: each dropped stim in droplist, d2: coord
                                    %put these into coords_outsample_alignment
                                    for idrop=1:ndrop
                                        istim=droplist{irept,iset}(idrop);
                                        coords_outsample_align{istim}(end+1,:)=coords_outsample_align_alldrop(idrop,:);
                                    end                                   
                                elseif ndrop>opts_xv_used.omit_per_fold %if all of the trial is dropped (corresponds to omit_per_fold=0)
                                    %transform entire repeat into consensus_insample as a reference, first doing a private pca
                                    rpca_outsample=reshape(resps(:,irept,:),[nstims,nrois_avail(subsamp_sel(iset))]);
                                    nonans_pca=find(all(~isnan(rpca_outsample),2));
                                    [u_nonan,s,v]=svd(rpca_outsample(nonans_pca,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
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
                                        coords_outsample_align{istim}(end+1,:)=coords_outsample_xform(istim,:);
                                    end
                                end % ndrop
                            end %irept
                        end %iset
                        %
                        % add to confusion matrices                
                        % confusion_matrices=zeros(nstims,nstims,ndecs,dmax,nsubs,npreprocs,nsubsamps_use); % d1: actual stim, d2: decoded stim, d3: decision rule, d4: dmax, d5: sub mean d6: normalize, d7: subsample set
                        for idec=1:ndecs
                            [confmtx,aux]=xval_confmtx_make(coords_outsample_align,coords_insample_align,dec_methods{idec});
                            dec_labels{idec}=aux.desc_brief;
                            confusion_matrices(:,:,idec,id,ipreproc,isub,isubsamp)=confusion_matrices(:,:,idec,id,ipreproc,isub,isubsamp)+confmtx;
                        end
                    end %ifold
                end %ixv_make
                if (if_3dplot) & (id==3)
                    hlid_rastim_trial_decode_3dplot;
                end
            end %dmax
        end %ipreproc
    end %isub (mean subtract)   
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
results.if_frozen=if_frozen;
results.if_debug=if_debug;
%
results.metadata=metadata;
results.dsids=dsids;
results.stims_avail=stims_avail; %list of available stimuli in each file, beginning at 1
results.rois_avail=rois_avail;
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
results.confusion_matrices_dims={'d1: actual stim, d2: decoded stim, d3: decision rule, d4: dmax, d5: sub mean d6: normalize, d7: subsample set'};
%
disp('results structure created.');
%
disp('confusion matrices can be plotted with hlid_rastim_trial_decode_confmtx')
