%hlid_mds_transform_jackstats analyzes the transformation between two representational spaces, 
% by jackknifing stimuli and datasets
% 
% Builds on hlid_geom_transform_stats, with plans for the following:
%  Makes use of nonstandard embeddings created by hlid_rastim_mds_coords_make
%    * consensuses across files are already calculated in hlid_rastim_mds_coords_make
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
%  Coordinates for all embeddings, for tntlabel and tnt3c, are first accumulated into z.[ref|adj]
%     z.[ref|adj].r_all.jackknife_omit_none{1} are the  non-jackknifed quantities
%     z.[ref|adj].r_all.jackknife_by_stim{:} have one stimulus omitted
%     z.[ref|adj].r_all.jackknife_by_file{:} have one set omitted
%
%      Consensuses are in  z.ref.r_all.jackknife_*{*}.data_knit{imeth,ism}.ds{1}
%
%
% makes use of code from
%  psg_majaxes, hlid_majaxes: examine axes identified in a transformation
%  psg_align_vara_demo: variance analysis after aligning multiple datasets grouped by condition
%  psg_geomodels_run: to determine transformation between two datasets
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
%  HLID_RASTIM_MDS_COORDS_DEMO, HLID_MDS_COORDS_GEOMODELS,
%  HLID_RASTIM_MDS_COORDS_MAKE, HLID_MDS_TRANSFORM_STATS.
%
hlid_setup;  %invoke hlid_localopts; set up opts_read and opts_plot
if ~exist('eig_tol') eig_tol=10^-6;end %tolerance for a negative eigenvalue
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
%get all coordinates
for ira=1:nra
    for imeth=1:nmeths
        for submean=0:if_submean
            %this is the consensus coordinate set, typically not in principal components
            coord_sets_all=z.(gp_labels{ira}).r_all.jackknife_omit_none{1}.data_knit{imeth,1+submean}.ds{1};
        end
    end
end %ira
%find max Euclidean dimension: survey the eigenvalues of the individual datasets to see where they first go negative
%
max_euc_dim=zeros(nfiles,nmeths,1+if_submean);
max_euc_dim_jackstim=zeros(nfiles,nmeths,1+if_submean,nstims);
for ifile=1:nfiles
    for imeth=1:nmeths
        for submean=0:if_submean
            if (ifile<=nfiles_ra(1))
                ira=1;
            else
                ira=2;
            end
            file_ptr=ifile-sum(nfiles_ra(1:ira-1));
            eivals=z.(gp_labels{ira}).r_all.jackknife_omit_none{1}.eivals{file_ptr,imeth,1+submean};
            max_euc_dim(ifile,imeth,1+submean)=-1+min(find(double([eivals;-1]>=-eig_tol)==0));
            for istim=1:nstims
                eivals_jackstim=z.(gp_labels{ira}).r_all.jackknife_by_stim{istim}.eivals{file_ptr,imeth,1+submean};
                max_euc_dim_jackstim(ifile,imeth,1+submean,istim)=-1+min(find(double([eivals_jackstim;-1]>=-eig_tol)==0));
            end
        end %if_submean
    end %imeth
end %ifile
%
%look at transformations between conditions, and look at magnif factor
%
%centering is taken care of by cosine and Pearson correlations:  Pearson is centered, cosine is not
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
ngps=2;
gps=[repmat(1,1,nfiles_ra(1)),repmat(2,1,nfiles_ra(2))];
gp_list={1:nfiles_ra(1),nfiles_ra(1)+[1:nfiles_ra(2)]};
nsets_gp=nfiles_ra;
nsets_gp_max=max(nsets_gp);
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
results.meths=meths;
results.meth_names_short=meth_names_short;
results.nmeths=nmeths;
results.meth_use_list=meth_use_list;
results.submean_use_list=submean_use_list;
%
results.nsubs=nsubs;
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
%
results.metadata=metadata;
results.dsids=dsids;
results.stimulus_names_display=stimulus_names_display;
%
results.nmodels=nmodels;
results.model_types_use=model_types_use;
results.opts_majaxes=opts_majaxes;
%
results.geo=cell(1+if_submean,nmeths,nembeds);
results.geo_majaxes=cell(1+if_submean,nmeths,nembeds);
results.geo_jack_by_stim=cell(1+if_submean,nmeths,nembeds);
results.geo_majaxes_jack_by_stim=cell(1+if_submean,nmeths,nembeds,nstims);
results.geo_majaxes_dims='d1: submean, d2: embedding method (mds,cosine,Pearson), d3: nembeds (1), d4: jack';
%
for imeth_ptr=1:length(meth_use_list)
    imeth=meth_use_list(imeth_ptr);
    for isubmean_ptr=1:length(submean_use_list)
        isubmean=submean_use_list(isubmean_ptr);
        for iembed=1:nembeds %do the embedding for each file, for this preprocessing choice
            %
            %fit geometric models listed in model_types_use
            %but no shuffling to determine singificance of nested models
            %
            for ijack=0:nstims
                d_ref=cell(1,max(dimlist));
                d_adj=cell(1,max(dimlist));
                for idim_ptr=1:length(dimlist)
                    npcs=dimlist(idim_ptr);
                    if ijack==0
                        d_ref{npcs}=z.ref.r_all.jackknife_omit_none{1}.data_knit{imeth,1+isubmean}.ds{1}{npcs};
                        d_adj{npcs}=z.adj.r_all.jackknife_omit_none{1}.data_knit{imeth,1+isubmean}.ds{1}{npcs};
                        stims_keep=[1:nstims];
                    else
                        d_ref{npcs}=z.ref.r_all.jackknife_by_stim{ijack}.data_knit{imeth,1+isubmean}.ds{1}{npcs};
                        d_adj{npcs}=z.adj.r_all.jackknife_by_stim{ijack}.data_knit{imeth,1+isubmean}.ds{1}{npcs};
                        stims_keep=setdiff(1:nstims,ijack);
                    end
                    if if_center
                        d_ref{npcs}=d_ref{npcs}-repmat(mean(d_ref{npcs},1),[size(d_ref{npcs},1),1]);
                        d_adj{npcs}=d_adj{npcs}-repmat(mean(d_adj{npcs},1),[size(d_adj{npcs},1),1]);
                    end
                end
                sa_ref=struct;
                sa_ref.typenames=stimulus_names_display(stims_keep);
                sa_adj=struct;
                sa_adj.typenames=stimulus_names_display(stims_keep);
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
                if ijack==0
                    results.geo{1+isubmean,imeth,iembed}=r_geo_all;
                    results.geo_majaxes{1+isubmean,imeth,iembed}=r_geo_majaxes;
                else
                    results.geo_jack_by_stim{1+isubmean,imeth,iembed,ijack}=r_geo_all;
                    results.geo_majaxes_jack_by_stim{1+isubmean,imeth,iembed,ijack}=r_geo_majaxes;
                end
            end %ijack
        end %iembed
    end %isubmean
end % imeth
%
%magnif factors, as a function of submean, method, dimension and jackknife
%
%d1: dim-1, d2: first then second magnif then lowest then geomean, d3: full then jackknife, d4: submean d5: method
magfacs=zeros(length(dimlist)-1,4,1+nstims,length(results.submean_use_list),length(results.meth_use_list));
%
for imeth_ptr=1:length(results.meth_use_list)
    imeth=meth_use_list(imeth_ptr);
    for isubmean_ptr=1:length(results.submean_use_list)
        isubmean=results.submean_use_list(isubmean_ptr);
        for ijack=0:nstims
            if (ijack==0)
                mf_all=results.geo_majaxes{1+isubmean,imeth};
            else
                mf_all=results.geo_majaxes_jack_by_stim{1+isubmean,imeth,1,ijack};
            end
            for k=2:dimlist(end)
                magfacs(k-1,1:2,ijack+1,1+isubmean,imeth)=mf_all{k,k}.ref.magnifs{1}(1:2)'; %take top two values
                magfacs(k-1,3,ijack+1,1+isubmean,imeth)=mf_all{k,k}.ref.magnifs{1}(end); %lowest value
                magfacs(k-1,4,ijack+1,1+isubmean,imeth)=geomean(mf_all{k,k}.ref.magnifs{1}); %lowest value
            end
        end
    end
end
%
%plots
%
figure;
set(gcf,'Position',[100 100 1200 800]);
set(gcf,'Name','maj axis ratios');
set(gcf,'NumberTitle','off');
for imeth_ptr=1:length(results.meth_use_list)
    imeth=meth_use_list(imeth_ptr);
    for isubmean_ptr=1:length(results.submean_use_list)
        isubmean=results.submean_use_list(isubmean_ptr);
        subplot(length(results.submean_use_list),length(results.meth_use_list),imeth_ptr+(isubmean_ptr-1)*length(results.meth_use_list));
        plot(dimlist(2:end),magfacs(:,1,1,1+isubmean,imeth)./magfacs(:,2,1,1+isubmean,imeth),'r'); %un-jackknifed ratio
        hold on;
        plot(dimlist(2:end),reshape(magfacs(:,1,2:end,1+isubmean,imeth)./magfacs(:,2,2:end,1+isubmean,imeth),length(dimlist)-1,nstims),'k'); %jackknifed ratio
        set(gca,'XLim',[1.5 results.dimlist(end)+0.5]);
        set(gca,'XTick',results.dimlist(2:end));
        xlabel('dim');
        set(gca,'YLim',[1 2]);
        set(gca,'YScale','log');
        ylabel('a1/a2');
        title(sprintf('sm=%1.0f %s',isubmean,results.meth_names_short{imeth}));
    end
end
disp('results structure created, consider saving it');
