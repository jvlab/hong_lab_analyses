%hlid_mds_transform_stats analyzes the transformation between two representational spaces, 
% including shuffles and resamplings for statistics
%
% Builds on hlid_geom_transform_stats, with plans for the following:
%  Makes use of nonstandard emvbeddings created by hlid_rastim_mds_coords_make
%    * specific for comparisons of tnt3c and tntlabel datasets
%      - assumes identical stimulus sets
%      - raw data read in hlid_rastim_mds_coords_make, and not here (so no call to hlid_rastim_trial_read)
%      - embedding calcs done in hlid_rastim_mds_coords_make
%    * trial-averaging is always done (via hlid_rastim_mds_coords_make)
%    * mean subtraction and normalization done as part of the nonstandard embeddings
%    * removed results fields not needed for further analysis
%    * since the data have already been read as a ref and adj group, the
%    order of the files in dsids and metadata is differnt c/w hlid_geom_transform_stats
%   
%  Jackknife via removing individual stimuli
%  Jackknife via removing individual preps
%
% makes use of code from
%  psg_majaxes, hlid_majaxes: examine axes identified in a transformation
%  psg_align_vara_demo: variance analysis after aligning multiple datasets grouped by condition
%  psg_geomodels_run: to determine transformation between two datasets
%
% Then shuffles the labels and re-analyzes; hlid_geom_transform_stats_summ
% will plot and summarize.  Main results saved in 'results'.
%
% The list of models to analyze is given by model_types. This should be a subset of
% the model_types field of psg_geomodels_define().  It defaults to 'affine_nooffset'
%
% If model_types is changed, then opts_majaxes.model_class_list should be
% changed to contain all of the model classes, typically {'affine','projective','pwaffine'}
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO,
%  PROCRUSTES_CONSENSUS, PROCRUSTES_COMPAT,
%  HLID_GEOM_TRANSFORM_STATS_3DPLOT, HLID_GEOM_TRANSFORM_STATS_SUMM, 
%  PSG_GEOMODELS_DEFINE, PSG_GEO_GENERAL,
%  HLID_RASTIM_TRIAL_DECODE, HLID_MAJAXES, PSG_ALIGN_VARA_DEMO, PSG_GEOMODELS_RUN, PSG_MAJAXES,
%  MULTI_SHUFF_GROUPS, MULTI_BOOT_GROUPS, HLID_GEOM_TRANSORM_STATS,
%  HLID_RASTIM_MDS_COORDS_DEMO, HLID_MDS_COORDS_GEOMODELS, HLID_RASTIM_MDS_COORDS_MAKE.
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
gp_labels={'ref','adj'}; %transformations from adjusted set into reference set, same as psg_geomodels_run
%
if ~exist('ref_coordfile') ref_coordfile='hlid_rastim_mds_make_kc-tnt3c_30Dec25.mat'; end
if ~exist('adj_coordfile') adj_coordfile='hlid_rastim_mds_make_kc-tntlabel_30Dec25.mat'; end
z=struct;
disp(sprintf(' loading ref coords from %s',ref_coordfile));
z.ref=load(ref_coordfile);
disp(z.ref.r_all);
disp(sprintf(' loading adj coords from %s',adj_coordfile));
z.adj=load(adj_coordfile);
disp(z.adj.r_all);
disp('all loaded')
disp(z);
%
nra=length(gp_labels);
dmax_avail=length(z.(gp_labels{1}).r_all.jackknife_omit_none{1}.data_knit{1,1}.ds{1});
%
nfiles_ra=zeros(1,nra);
for ira=1:nra
    nfiles_ra(ira)=length(z.(gp_labels{ira}).r_all.jackknife_omit_none{1}.filenames_short);
    disp(sprintf('files for %s: %2.0f',gp_labels{ira},nfiles_ra(ira)))
end
nfiles=sum(nfiles_ra);
disp(sprintf('number of files: %4.0f',nfiles));
%
stimulus_names_display=z.ref.r_all.jackknife_omit_none{1}.data{1,1}.sas{1}.typenames;
nstims=length(stimulus_names_display);
disp(sprintf('number of stimuli: %4.0f',nstims));
%
if_submean=1;
meths=z.ref.r_all.jackknife_omit_none{1}.meths;
nmeths=length(meths);
meth_names_short=cell(1,nmeths);
disp(sprintf('number of embedding methods: %4.0f',nmeths));
for imeth=1:nmeths
    meth_names_short{imeth}=meths{imeth}.name_short;
    disp(sprintf(' method %2.0f: %30s transform: %20s (brief: %20s or %12s)',imeth,meths{imeth}.name_full,meths{imeth}.xform,meths{imeth}.name_short,meths{imeth}.name_file));
end
meth_use_list=getinp('list to use','d',[1 nmeths],[1:nmeths]);
submean_use_list=getinp('list of subtract-mean options','d',[0 1],[0 1]);
%
dimpair_opts={'ref and adj have same dimensions','ref and adj step through all dimensions','ref is always <= adj','ref is always >= adj'};
if ~exist('if_econ_svd') if_econ_svd=1; end
%
% here is where we survey the data for a summary of frac variance explained by dimension, with error bars based on jackknifes
% compare participation ratio across jackknife sets and methods
% eigenvalues might be negative for the individual coordinate sets prior to
% forming consensus; should only look at particpation ratio and power ratio until some of
% the eigenvalues of the individual datasets become negative
%
for ira=1:nra
    for imeth=1:nmeths
        for submean=0:if_submean
            %this is the consensus coordinate set, typically not principal components. Need to do PCA to get eigenvalues and participation ratio
            %Since we are doing PCA, power ratio will always be 1.0.
            coord_sets_all=z.(gp_labels{ira}).r_all.jackknife_omit_none{1}.data_knit{imeth,1+submean}.ds{1};
        end
    end
end %ira
%look at transformations between conditions, and look at magnif factor
%stats
%
%default is to not allow offset and center data;
%alternately could allow for offset and/or not center
%
model_types_def=psg_geomodels_define();
model_types_use=struct;
if ~exist('model_types') model_types={'affine_nooffset'}; end
if_cycle=1; %only relevant for projective models
model_types_use.model_types=model_types;
nmodels=length(model_types);
for imodel=1:nmodels
    mname=model_types{imodel};
    model_types_use.(mname)=model_types_def.(mname);
    if contains(model_types_use.(mname).class,'projective')
        model_types_use.(mname).opts.if_cycle=if_cycle;
    end
end
%
if ~exist('if_center') if_center=1; end %whether to center the data after forming the consensus
%
if ~exist('opts_majaxes') opts_majaxes=struct; end
opts_majaxes=filldefault(opts_majaxes,'model_class_list',{'affine'}); %only affine models are calculated
%
%
if ~exist('opts_pcon') opts_pcon=struct; end %consensus options
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
%
%differences c/w hlid_geom_transform_stats
embed_labels={'embed trial-averages'}; % only embed trial averages
nembeds=length(embed_labels); %here, =1
sub_labels={''}; %mean subtraction taken care of at embedding stage
nsubs=length(sub_labels); %here, =1
preproc_labels={'raw'}; %taken care of at embedding stage
npreprocs=length(preproc_labels); %here,=1
%
% fields reconstructed from files already read
dsids=cell(nfiles,1);
dsids(1:nfiles_ra(1))=z.ref.r_all.jackknife_omit_none{1}.filenames_short;
dsids(nfiles_ra(1)+[1:nfiles_ra(2)])=z.adj.r_all.jackknife_omit_none{1}.filenames_short;
metadata=cell(nfiles,1);
for ifile=1:nfiles_ra(1)
    metadata{ifile}=z.ref.r_all.jackknife_omit_none{1}.data{1,1}.sets{ifile};
end
for ifile=1:nfiles_ra(2)
    metadata{nfiles_ra(1)+ifile}=z.adj.r_all.jackknife_omit_none{1}.data{1,1}.sets{ifile};
end
%
if ~exist('dmax') 
    if if_debug
        dmax=min(5,dmax_avail);
        dimlist_def=[1:4];
    else
        dmax=min(nstims-1,dmax_avail); %max representational space to create
        dimlist_def=[2:6];
    end
end
if ~exist('dmin') 
    dmin=1;
end
if_ok=0;
while (if_ok==0)
    dimlist=getinp('list of dimensions for adj and ref','d',[dmin dmax],dimlist_def);
    for k=1:length(dimpair_opts)
        disp(sprintf('%1.0f->%s',k,dimpair_opts{k}));
    end
    dimpair_opt=getinp('choice','d',[1 length(dimpair_opts)]);
    switch dimpair_opt
        case 1
            dimpair_list=repmat(dimlist(:),1,2);
        case {2,3,4}
            dimpair_list=zeros(0,2)
            for iref=1:length(dimlist)
                for iadj=1:length(dimlist)
                    if dimpair_opt==2 | (dimpair_opt==3 & dimlist(iref)<=dimlist(iadj)) | (dimpair_opt==4 & dimlist(iref)>=dimlist(iadj))
                        dimpair_list(end+1,:)=[dimlist(iref) dimlist(iadj)];
                    end
                end
            end
    end
    disp('dimensions to be examined');
    for igp=1:length(gp_labels)
        disp(sprintf('%4s, %s',gp_labels{igp},sprintf('%3.0f ',dimpair_list(:,igp))));
    end
    if_ok=getinp('1 if ok','d',[0 1]);
end
%
%get grouping information: group 1 is ref, group 2 is adj
%
nsets=nfiles;
sets=cell(1,nsets);
for iset=1:nsets
    sets{iset}.label=dsids{iset};
    disp(sprintf(' file (set) %3.0f: %s',iset,sets{iset}.label));
end
for igp=1:length(gp_labels)
    disp(sprintf(' group %1.0f is %s',igp,gp_labels{igp}))
end
[ngps,gps,gp_list,nsets_gp,nsets_gp_max]=psg_getgps(sets,length(gp_labels));
%
opts_pcon.allow_scale=getinp('1 to allow scaling within groups','d',[0 1],0);
scaling_token=(opts_pcon.allow_scale==1);
if opts_pcon.allow_reflection==1
    reflection_token='best';
else
    reflection_token=false;
end
%
if_3dplot=getinp('1 to plot response space for 3d consensus','d',[0 1]);
%
%get between-group shuffle information
%
opts_shuff=struct;
opts_shuff.if_ask=-1;
opts_shuff.if_reduce=0;
if if_debug
    opts_shuff.nshuffs=10;
else
    opts_shuff.nshuffs=1000;
end
[shuffs_between,gp_info,opts_shuff_used]=multi_shuff_groups(gps,opts_shuff);
nshuffs_between=size(shuffs_between,1);
shuff_gp_selects=cell(1,ngps);
shuff_gp_origs=cell(1,ngps);
for igp=1:ngps
    shuff_gp_selects{igp}=zeros(nshuffs_between,nsets_gp(igp)); %which sets are in each group
    shuff_gp_origs{igp}=zeros(nshuffs_between,nsets_gp(igp)); %original group membership in each group
end
if_keep_all_shuff=0;
if nshuffs_between>0
    if_keep_all_shuff=getinp('1 to keep all outputs from shuffles, 0 for outputs only needed for stats','d',[0 1],if_debug);
    for ishuff=1:nshuffs_between
        gp_select=cell(1,ngps);
        gp_orig=cell(1,ngps);
        for igp=1:ngps
            gp_select{igp}=find(gps(shuffs_between(ishuff,:))==igp); %shuffled datasets for group igp
            gp_orig{igp}=gps(gp_select{igp});
            shuff_gp_selects{igp}(ishuff,:)=gp_select{igp};
            shuff_gp_origs{igp}(ishuff,:)=gp_orig{igp};
        end
    end
end
%
% get bootstrap settings
%
opts_boot=opts_shuff;
[boots_within,gp_info_boot,opts_boot_used]=multi_boot_groups(gps,opts_boot);
nboots_within=size(boots_within,1);
boot_gp_selects=cell(1,ngps);
boot_gp_origs=cell(1,ngps);
for igp=1:ngps
    boot_gp_selects{igp}=zeros(nboots_within,nsets_gp(igp)); %which sets are in each group
    boot_gp_origs{igp}=zeros(nboots_within,nsets_gp(igp)); %original group membership in each group
end
if nboots_within>0
    if_keep_all_boot=getinp('1 to keep all outputs from bootstraps, 0 for outputs only neeed for stats','d',[0 1],if_debug);
    for iboot=1:nboots_within
        gp_select=cell(1,ngps);
        gp_orig=cell(1,ngps);
        for igp=1:ngps
            gp_select{igp}=boots_within(iboot,gp_list{igp});
            gp_orig{igp}=gps(gp_select{igp});
            boot_gp_selects{igp}(iboot,:)=gp_select{igp};
            boot_gp_origs{igp}(iboot,:)=gp_orig{igp}; %should always = igp
        end
    end
end
for igp=1:ngps
    if any(boot_gp_origs{igp}~=igp)
        warning(sprintf(' bootstrap cross-check fails for group %1.0f',igp));
    end
end
%
%save settings
%
results=struct;
results.nstims=nstims;
results.nfiles=nfiles;
results.nsets=nsets;
%
results.gp_labels=gp_labels;
results.ngps=ngps;
results.gps=gps;
results.gp_list=gp_list;
results.nsets_gp=nsets_gp;
%
results.nsubs=nsubs;
results.meths=meths;
results.meth_names_short=meth_names_short;
results.nmeths=nmeths;
results.meth_use_list=meth_use_list;
%
results.sub_labels=sub_labels;
results.npreprocs=npreprocs;
results.preproc_labels=preproc_labels;
results.nembeds=nembeds;
results.embed_labels=embed_labels;
results.dmax=dmax;
results.dmin=dmin;
results.dimlist=dimlist;
results.dimpair_opts=dimpair_opts;
results.dimpair_opt=dimpair_opt;
results.dimpair_list=dimpair_list;
results.if_frozen=if_frozen;
results.if_debug=if_debug;
results.if_econ_svd=if_econ_svd;
%
results.metadata=metadata;
results.dsids=dsids;
results.stimulus_names_display=stimulus_names_display;
%
results.nmodels=nmodels;
results.model_types_use=model_types_use;
results.opts_majaxes=opts_majaxes;
%
results.opts_pcon=opts_pcon;
%
results.nshuffs_between=nshuffs_between;
results.shuffs_between=shuffs_between;
results.shuff_gp_selects=shuff_gp_selects;
results.shuff_gp_origs=shuff_gp_origs;
results.opts_shuff=opts_shuff_used;
%
results.nboots_within=nboots_within;
results.boots_within=boots_within;
results.boot_gp_selects=boot_gp_selects;
results.opts_boot=opts_boot_used;
%
consensus_init=cell(1,length(dimlist)); %initial guess, uniform across nsubs, npreprocs, nembeds, if_submean, nmeths
%
results.geo=cell(1+if_submean,nmeths,nembeds);
results.geo_majaxes=cell(1+if_submean,nmeths,nembeds);
results.geo_majaxes_dims='d1: submean, d2: embedding method (mds,cosine,Pearson), d3: nembeds (1), d4: nshuffs';
if if_keep_all_shuff
    results.geo_shuff=cell(1+if_submean,nmeths,nembeds,nshuffs_between);
end
results.geo_majaxes_shuff=cell(1+if_submean,nmeths,nembeds,nshuffs_between);
%
for imeth_ptr=1:length(meth_use_list)
    imeth=meth_use_list(imeth_ptr);
    for isubmean_ptr=1:length(submean_use_list)
        isubmean=submean_use_list(isubmean_ptr);
        coords_embedded=cell(nembeds,length(dimlist));
        consensus_glbl=cell(nembeds,length(dimlist)); %consensus across all groups
        consensus_bygp=cell(nembeds,length(dimlist),ngps); %consensus within groups
        znew_glbl=cell(nembeds,length(dimlist)); %coords after alignment to global consensus
        znew_bygp=cell(nembeds,length(dimlist)); %coords after alignment to in-group consensus
        for iembed=1:nembeds %do the embedding for each file, for this preprocessing choice
            npts_embed=nstims;
            for idim_ptr=1:length(dimlist)
                npcs=dimlist(idim_ptr);
                coords_embedded{iembed,idim_ptr}=zeros(npts_embed,npcs,nfiles);
                %retrieve coordinates from r_all, calculated in hlid_mds_coords_make
                for ifile=1:nfiles
                    ira=gps(ifile);
                    ifile_ptr=ifile-sum(nfiles_ra(1:ira-1));
                    coords_embedded{iembed,idim_ptr}(:,:,ifile)=z.(gp_labels{ira}).r_all.jackknife_omit_none{1}.data{imeth,1+isubmean}.ds{ifile_ptr}{idim_ptr};
                end
                %
                %consensus within groups, without shuffling, and use global consensus for iembed=1 for initial guess and alignment once it is created
                %
                if isempty(consensus_init{idim_ptr})
                    [consensus_glbl{iembed,idim_ptr},znew_glbl{iembed,idim_ptr}]=procrustes_consensus(coords_embedded{iembed,idim_ptr},opts_pcon);
                    consensus_init{1,idim_ptr}=consensus_glbl{iembed,idim_ptr};
                else
                    opts_pcon_use.initialize_set=0;
                    opts_pcon_use.initial_guess=consensus_init{1,idim_ptr};
                    [consensus_glbl{iembed,idim_ptr},znew_glbl{iembed,idim_ptr}]=procrustes_consensus(coords_embedded{iembed,idim_ptr},opts_pcon_use);
                end
                disp(sprintf(' consensus made across groups,      dim %2.0f, isubmean=%1.0f, method %s',npcs,isubmean,meth_names_short{imeth}));
                znew_bygp{iembed,idim_ptr}=zeros(npts_embed,npcs,nfiles); %fill third coord according to which group
                for igp=1:ngps
                    opts_pcon_use=opts_pcon;
                    opts_pcon_use.initialize_set=0;
                    opts_pcon_use.initial_guess=consensus_init{1,idim_ptr};
                    [consensus_bygp{iembed,idim_ptr,igp},znew_bygp{iembed,idim_ptr}(:,:,gp_list{igp})]=procrustes_consensus(coords_embedded{iembed,idim_ptr}(:,:,gp_list{igp}),opts_pcon_use);
                    disp(sprintf(' consensus made for group %1.0f (%4s), dim %2.0f, isubmean=%1.0f, method %s',igp,gp_labels{igp},npcs,isubmean,meth_names_short{imeth})); 
                end
                if (npcs==3) & if_3dplot
                    title_string=sprintf('isubmean=%1.0f, method %s',isubmean,meth_names_short{imeth});
                    hlid_geom_transform_stats_3dplot;
                end
            end %idim_ptr
            %
            %fit geometric models listed in model_types_use
            %but no shuffling to determine singificance of nested models
            %
            d_ref=cell(1,max(dimlist));
            d_adj=cell(1,max(dimlist));
            for idim_ptr=1:length(dimlist)
                npcs=dimlist(idim_ptr);
                d_ref{npcs}=consensus_bygp{iembed,idim_ptr,1}; %unshuffled
                d_adj{npcs}=consensus_bygp{iembed,idim_ptr,2}; %unshuffled
                if if_center
                    d_ref{npcs}=d_ref{npcs}-repmat(mean(d_ref{npcs},1),[npts_embed,1]);
                    d_adj{npcs}=d_adj{npcs}-repmat(mean(d_adj{npcs},1),[npts_embed,1]);
                end
            end
            sa_ref=struct;
            sa_ref.typenames=stimulus_names_display;
            sa_adj=struct;
            sa_adj.typenames=stimulus_names_display;
            %
            r_geo_all=cell(max(dimlist),max(dimlist)); %in keeping with psg_geomodels_run
            r_geo_base=cell(max(dimlist),max(dimlist)); %for use for shuffles
            r_geo=struct;
            for idim_pair=1:size(dimpair_list,1)
                ref_dim=dimpair_list(idim_pair,1);
                adj_dim=dimpair_list(idim_pair,2);
                r_geo.model_types_def=model_types_use;
                r_geo.ref_dim=ref_dim;
                r_geo.adj_dim=adj_dim;
                r_geo.ref_file=cat(2,sprintf('consensus of %2.0f files:,',length(gp_list{1})),sets{gp_list{1}(1)}.label,'...',sets{gp_list{1}(end)}.label);
                r_geo.adj_file=cat(2,sprintf('consensus of %2.0f files:,',length(gp_list{2})),sets{gp_list{2}(1)}.label,'...',sets{gp_list{2}(end)}.label);
                r_geo.if_center=if_center;
                r_geo.if_cycle=if_cycle;
                r_geo.nestdim_list=[];
                r_geo.d=NaN(nmodels,1);
                r_geo.adj_model=cell(nmodels,1); %not in psg_geo_models_run but could be useful (model fit)
                r_geo.transforms=cell(nmodels,1);
                r_geo.opts_model_used=cell(nmodels,1);
                r_geo_base{ref_dim,adj_dim}=r_geo;
                for imodel=1:nmodels
                    model_type=model_types{imodel};
                    model_desc_string=sprintf('model type: %s, adj dim %2.0f ref dim %2.0f',model_type,adj_dim,ref_dim);
                    disp(model_desc_string);
                    if adj_dim<model_types_def.(model_type).min_inputdims
                        disp(sprintf('model type %s skipped, requires input dimension of at least %2.0f',model_type,model_types_def.(model_type).min_inputdims));
                    else
                        opts_model=model_types_def.(model_type).opts;
                        model_class=model_types_def.(model_type).class;
                        [r_geo.d(imodel),r_geo.adj_model{imodel},r_geo.transforms{imodel},r_geo.opts_model_used{imodel}]=...
                            psg_geo_general(d_ref{ref_dim},d_adj{adj_dim},model_class,opts_model);
                    end
                end %imodel
                r_geo_all{ref_dim,adj_dim}=r_geo;
            end %idim_pair
            [r_geo_majaxes,opts_majaxes_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,r_geo_all,opts_majaxes);
            results.geo{1+if_submean,imeth,iembed}=r_geo_all;
            results.geo_majaxes{1+if_submean,imeth,iembed}=r_geo_majaxes;
            %
            %shuffles
            %
            if nshuffs_between>0
                for ishuff=1:nshuffs_between
                    r_geo_all_shuff=cell(max(dimlist),max(dimlist)); %in keeping with psg_geomodels_run
                    d_ref_shuff=cell(1,max(dimlist));
                    d_adj_shuff=cell(1,max(dimlist));
                    %create consensus within groups, after shuffling between groups, using global consensus for initialization and alignment
                    opts_pcon_use=opts_pcon;
                    opts_pcon_use.initialize_set=0;
                    consensus_shuff_bygp=cell(length(dimlist),ngps); %compute the consensus for each group, each dimension, this embedding, and center if necessary
                    for idim_ptr=1:length(dimlist)
                        for igp=1:ngps
                            opts_pcon_use.initial_guess=consensus_init{1,idim_ptr};
                            consensus_shuff_bygp{idim_ptr,igp}=procrustes_consensus(coords_embedded{iembed,idim_ptr}(:,:,shuff_gp_selects{igp}(ishuff,:)),opts_pcon_use);
                            if if_center
                                consensus_shuff_bygp{idim_ptr,igp}=consensus_shuff_bygp{idim_ptr,igp}-repmat(mean(consensus_shuff_bygp{idim_ptr,igp},1),[npts_embed,1]);
                            end
                            d_ref_shuff{dimlist(idim_ptr)}=consensus_shuff_bygp{idim_ptr,1};
                            d_adj_shuff{dimlist(idim_ptr)}=consensus_shuff_bygp{idim_ptr,2};
                        end %igp
                    end %idim_ptr
                    for idim_pair=1:size(dimpair_list,1)
                        ref_dim=dimpair_list(idim_pair,1);
                        adj_dim=dimpair_list(idim_pair,2);
                        r_geo_shuff=r_geo_base{ref_dim,adj_dim};
                        r_geo_shuff.ref_file=sprintf('ref, shuff_between %5.0f of %5.0f',ishuff,nshuffs_between);
                        r_geo_shuff.adj_file=sprintf('adj, shuff_between %5.0f of %5.0f',ishuff,nshuffs_between);
                        idim_ref_ptr=find(dimlist==ref_dim);
                        idim_adj_ptr=find(dimlist==adj_dim);
                        %
                        r_geo_shuff.d=NaN(nmodels,1);
                        r_geo_shuff.adj_model=cell(nmodels,1); %not in psg_geo_models_run but could be useful (model fit)
                        r_geo_shuff.transforms=cell(nmodels,1);
                        r_geo_shuff.opts_model_used=cell(nmodels,1);
                        %fit the model
                        for imodel=1:nmodels
                            model_type=model_types{imodel};
                            if adj_dim>=model_types_def.(model_type).min_inputdims
                                opts_model=model_types_def.(model_type).opts;
                                model_class=model_types_def.(model_type).class;
                                [r_geo_shuff.d(imodel),r_geo_shuff.adj_model{imodel},r_geo_shuff.transforms{imodel},r_geo_shuff.opts_model_used{imodel}]=...
                                    psg_geo_general(d_ref_shuff{ref_dim},d_adj_shuff{adj_dim},model_class,opts_model);
                            end
                        end %imodel
                         r_geo_all_shuff{ref_dim,adj_dim}=r_geo_shuff;
                    end %idim_pair
                    %determine major axes
                    [r_geo_majaxes_shuff,opts_shuff_majaxes_used]=psg_majaxes(d_ref_shuff,sa_ref,d_adj_shuff,sa_adj,r_geo_all_shuff,opts_majaxes);
                    if if_keep_all_shuff
                        results.geo_shuff{1+if_submean,imeth,iembed,ishuff}=r_geo_all_shuff;
                    else %strip fields from r_geo_majaxes_shuff
                        r_geo_majaxes_shuff_full=r_geo_majaxes_shuff;
                        r_geo_majaxes_shuff=cell(size(r_geo_majaxes_shuff_full));
                        for ird=1:size(r_geo_majaxes_shuff,1)
                            for iad=1:size(r_geo_majaxes_shuff,2);
                                if ~isempty(r_geo_majaxes_shuff_full{ird,iad})
                                    r_geo_majaxes_shuff{ird,iad}.ref.magnifs=r_geo_majaxes_shuff_full{ird,iad}.ref.magnifs;
                                    r_geo_majaxes_shuff{ird,iad}.ref.magnif_ratio=r_geo_majaxes_shuff_full{ird,iad}.ref.magnif_ratio;
                                end %empty
                            end %iad
                        end %ird
                    end
                    results.geo_majaxes_shuff{1+if_submean,imeth,ishuff}=r_geo_majaxes_shuff;
                end %ishuff
                disp(sprintf(' %5.0f shuffles  between ref and adj done',nshuffs_between));
            end %nshuffs_between
            %
            %bootstraps
            %
            if nboots_within>0
                for iboot=1:nboots_within
                    r_geo_all_boot=cell(max(dimlist),max(dimlist)); %in keeping with psg_geomodels_run
                    d_ref_boot=cell(1,max(dimlist));
                    d_adj_boot=cell(1,max(dimlist));
                    %create consensus within groups, after bootstraps within groups, using global consensus for initialization and alignment
                    opts_pcon_use=opts_pcon;
                    opts_pcon_use.initialize_set=0;
                    consensus_boot_bygp=cell(length(dimlist),ngps); %compute the consensus for each group, each dimension, this embedding, and center if necessary
                    for idim_ptr=1:length(dimlist)
                        for igp=1:ngps
                            opts_pcon_use.initial_guess=consensus_init{1,idim_ptr};
                            consensus_boot_bygp{idim_ptr,igp}=procrustes_consensus(coords_embedded{iembed,idim_ptr}(:,:,boot_gp_selects{igp}(iboot,:)),opts_pcon_use);
                            if if_center
                                consensus_boot_bygp{idim_ptr,igp}=consensus_boot_bygp{idim_ptr,igp}-repmat(mean(consensus_boot_bygp{idim_ptr,igp},1),[npts_embed,1]);
                            end
                            d_ref_boot{dimlist(idim_ptr)}=consensus_boot_bygp{idim_ptr,1};
                            d_adj_boot{dimlist(idim_ptr)}=consensus_boot_bygp{idim_ptr,2};
                        end %igp
                    end %idim_ptr
                    for idim_pair=1:size(dimpair_list,1)
                        ref_dim=dimpair_list(idim_pair,1);
                        adj_dim=dimpair_list(idim_pair,2);
                        r_geo_boot=r_geo_base{ref_dim,adj_dim};
                        r_geo_boot.ref_file=sprintf('ref, boot_within %5.0f of %5.0f',iboot,nboots_within);
                        r_geo_boot.adj_file=sprintf('adj, boot_within %5.0f of %5.0f',iboot,nboots_within);
                        idim_ref_ptr=find(dimlist==ref_dim);
                        idim_adj_ptr=find(dimlist==adj_dim);
                        %
                        r_geo_boot.d=NaN(nmodels,1);
                        r_geo_boot.adj_model=cell(nmodels,1); %not in psg_geo_models_run but could be useful (model fit)
                        r_geo_boot.transforms=cell(nmodels,1);
                        r_geo_boot.opts_model_used=cell(nmodels,1);
                        %fit the model
                        for imodel=1:nmodels
                            model_type=model_types{imodel};
                            if adj_dim>=model_types_def.(model_type).min_inputdims
                                opts_model=model_types_def.(model_type).opts;
                                model_class=model_types_def.(model_type).class;
                                [r_geo_boot.d(imodel),r_geo_boot.adj_model{imodel},r_geo_boot.transforms{imodel},r_geo_boot.opts_model_used{imodel}]=...
                                    psg_geo_general(d_ref_boot{ref_dim},d_adj_boot{adj_dim},model_class,opts_model);
                            end
                        end %imodel
                         r_geo_all_boot{ref_dim,adj_dim}=r_geo_boot;
                    end %idim_pair
                    %determine major axes
                    [r_geo_majaxes_boot,opts_boot_majaxes_used]=psg_majaxes(d_ref_boot,sa_ref,d_adj_boot,sa_adj,r_geo_all_boot,opts_majaxes);
                    if if_keep_all_boot
                        results.geo_boot{1+if_submean,imeth,iembed,iboot}=r_geo_all_boot;
                    else %strip fields from r_geo_majaxes_boot
                        r_geo_majaxes_boot_full=r_geo_majaxes_boot;
                        r_geo_majaxes_boot=cell(size(r_geo_majaxes_boot_full));
                        for ird=1:size(r_geo_majaxes_boot,1)
                            for iad=1:size(r_geo_majaxes_boot,2);
                                if ~isempty(r_geo_majaxes_boot_full{ird,iad})
                                    r_geo_majaxes_boot{ird,iad}.ref.magnifs=r_geo_majaxes_boot_full{ird,iad}.ref.magnifs;
                                    r_geo_majaxes_boot{ird,iad}.ref.magnif_ratio=r_geo_majaxes_boot_full{ird,iad}.ref.magnif_ratio;
                                    r_geo_majaxes_boot{ird,iad}.ref.projections=r_geo_majaxes_boot_full{ird,iad}.ref.projections; %also keep projections, for confidence limits
                                    %
                                    r_geo_majaxes_boot{ird,iad}.adj.magnifs=r_geo_majaxes_boot_full{ird,iad}.adj.magnifs;
                                    r_geo_majaxes_boot{ird,iad}.adj.magnif_ratio=r_geo_majaxes_boot_full{ird,iad}.adj.magnif_ratio;
                                    r_geo_majaxes_boot{ird,iad}.adj.projections=r_geo_majaxes_boot_full{ird,iad}.adj.projections; %also keep projections, for confidence limits
                                end %empty
                            end %iad
                        end %ird
                    end
                    results.geo_majaxes_boot{1+if_submean,imeth,iembed,iboot}=r_geo_majaxes_boot;
                end %iboot
                disp(sprintf(' %5.0f bootstraps within ref and adj done',nboots_within));
            end %nboots_within
        end %iembed
    end %isubmean
end % imeth
%
disp('results structure created, consider saving it');
disp('summary plots can be made with hlid_geom_transform_stats_summ')
