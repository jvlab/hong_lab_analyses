%hlid_majaxes_eval2: further anlaysis of transformations of an affine geometric model to determine major axes
% for Hong Lab data
%
% derived from hlid_majaxes, designed to look at transformation between ORN space and KC space, 
%   compare the major axes of the transformation with the PC's of the space
%   compare the sqrt(variance) on each PC or principal axis with the magnif factors 
%
% Note: This does not check that stimuli are in the same order as used for the geometric modeling in psg_geomodels_run,
% which identifies the stimuli in common and then analyzes them in alphabetical order.
% For most purposes, the stimulus order needs to be set as follows:
% hlid_setup
% opts_read, opts_plot, opts_multim_def initialized for hlid data
% opts_majaxes.plot_order=display_orders.kcmerge
% opts_majaxes = 
%  struct with fields:
%    plot_order: {'2h'  'IaA'  'pa'  '2-but'  'eb'  'ep'  '6al'  't2h'  '1-8ol'  '1-5ol'  '1-6ol'  'PAA'  'ms'  'B-myr'  'euc'  '-aPine'  'pfo'}
%    plot_pairs: []
%
% compared to hlid_majaxes_eval:
%    *pccs and pcen don't have an ipw argument (but still are recomputed)
%    *computes stats on modeled transformation of adjusted dataset
%    *Some outputs reformatted
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE, HLID_SETUP, PSG_MAJAXES, HLID_MAJAXES_COMPARE,
%   PSG_MAJAXES_REORDER, PSG_PCA_OFFSET, PSG_GEOMODELS_APPLY, PSG_GEOMODELS_DEFINE.
%
% 30Jun25: fix bug: only consider stimuli in common to both ref and adj in computing statistics
%
hlid_setup;
model_types_def=psg_geomodels_define;
%
if ~exist('adj_ref_dims') 
    adj_ref_dims=[2 3 4 5];
end
if ~exist('adj_ref_pairs')
    adj_ref_pairs=repmat(adj_ref_dims',1,2);
end
if ~exist('tol')
    tol=10^-5;
end
%
if ~exist('opts_read') opts_read=struct;end
opts_read=filldefault(opts_read,'if_log',1);
opts_read=filldefault(opts_read,'if_justsetup',0);
opts_read=filldefault(opts_read,'data_fullname_def',hlid_opts.coord_data_fullname_def);
opts_read=filldefault(opts_read,'setup_fullname_def',[]);
opts_read=filldefault(opts_read,'domain_list',hlid_opts.domain_list_def);
opts_read=filldefault(opts_read,'setup_suffix',hlid_opts.setup_setup_suffix);
opts_read=filldefault(opts_read,'type_class_aux',hlid_opts.type_class_def);
%
if ~exist('opts_majaxes') opts_majaxes=struct;end
opts_majaxes=filldefault(opts_majaxes,'plot_pairs',[]);
opts_majaxes=filldefault(opts_majaxes,'plot_order',display_orders.kcclust);
%
fn_geo=getinp('file name containing results of psg_geomodels_run','s',[]);
results_geo=getfield(load(fn_geo),'results');
%
fn_ref=[];
fn_adj=[];
for k=1:size(results_geo,1)
    for m=1:size(results_geo,2)
        if ~isempty(results_geo{k,m})
            if isempty(fn_ref)
                fn_ref=results_geo{k,m}.ref_file;
            end
            if isempty(fn_adj)
                fn_adj=results_geo{k,m}.adj_file;
            end
        end
    end
end
fn_ref=getinp('file name for dataset that is the reference','s',[],fn_ref);
[d_ref,sa_ref]=psg_read_coorddata(fn_ref,[],opts_read);
%
fn_adj=getinp('file name for dataset that is adjusted','s',[],fn_adj);
[d_adj,sa_adj]=psg_read_coorddata(fn_adj,[],opts_read);
%
[results_axes,opts_majaxes_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts_majaxes);
%
%find overlap set, and remove non-overlaps
%
ref_ovlps=zeros(sa_ref.nstims,1);
for istim=1:sa_ref.nstims
    ref_ovlps(istim)=double(length(strmatch(sa_ref.typenames{istim},sa_adj.typenames,'exact'))==1);
end
disp(sprintf('ref set has %3.0f stimuli, %3.0f in common with adj',sa_ref.nstims,sum(ref_ovlps)));
adj_ovlps=zeros(sa_adj.nstims,1);
for istim=1:sa_adj.nstims
    adj_ovlps(istim)=double(length(strmatch(sa_adj.typenames{istim},sa_ref.typenames,'exact'))==1);
end
disp(sprintf('adj set has %3.0f stimuli, %3.0f in common with ref',sa_adj.nstims,sum(adj_ovlps)));
%
typenames_ovlps=sa_ref.typenames(ref_ovlps>0); %typenames that overlap
%
stim_select=cell(0);
plot_order=opts_majaxes_used.plot_order;
for istim=1:length(plot_order)
    if strmatch(plot_order{istim},sa_adj.typenames,'exact')>0 & strmatch(plot_order{istim},sa_ref.typenames,'exact')>0
        stim_select{end+1}=plot_order{istim};
    end
end
%
adj_ref_labels={'adj','ref'};
nar=length(adj_ref_labels);
%
%extract magnifications and vectors, and reorder/subtract mean to allow direct comparison to plots in psg_majaxes
%
n_pairs=size(adj_ref_pairs,1);
magnifs=cell(n_pairs,1); %magnifications
pccs=cell(n_pairs,1); %principal component coordinates of data used in psg_geomodels_run
pcen=cell(n_pairs,1); %principal components recomputed after centering
prjs=cell(n_pairs,1); %principal eigenvectors of the transformation, projected into odor space
pamc=cell(n_pairs,1); %principal comopnents of modeled transformaton of adj dataset, after centering
for iap=1:n_pairs
    id_adj=adj_ref_pairs(iap,1);
    id_ref=adj_ref_pairs(iap,2);
    if ~isempty(results_axes{id_ref,id_adj})
        ref_dim=results_axes{id_ref,id_adj}.ref_dim;
        adj_dim=results_axes{id_ref,id_adj}.adj_dim;
        disp(' ');
        disp(sprintf(' adj dim %2.0f ref dim %2.0f',adj_dim,ref_dim));
        magnifs{iap}=cell(0); %magnifications
        pccs{iap}=cell(0); %pcs of data
        pcen{iap}=cell(0); %pcs recomputed after recentering
        prjs{iap}=cell(0); %principal eigenvectors of the transformation, projected into odor space
        %
        for im_ptr=1:length(results_axes{id_ref,id_adj}.model_types) %step through each model
            model_type=results_axes{id_ref,id_adj}.model_types{im_ptr};
            npw=size(results_axes{id_ref,id_adj}.adj.eivecs{im_ptr},3); %number of transformations
            %awkward code; pccs, pcen do not depend on ipw
            for ipw=1:npw %step through each transformation
                disp(sprintf(' model type %20s  component %1.0f',model_type,ipw));
                magnifs{iap}{im_ptr,ipw}=struct;
                pccs{iap}{im_ptr}=struct; %coordinates from data, as used by psg_geomodels_run:  pc's computed around zero, but then centered
                pcen{iap}{im_ptr}=struct; %pcs after centering
                prjs{iap}{im_ptr,ipw}=struct; %projections on princiapl axes of transformation
                %stdvs around zero
                pccs_stdexp0{iap}{im_ptr}=struct;
                pcen_stdexp0{iap}{im_ptr}=struct;
                prjs_stdexp0{iap}{im_ptr,ipw}=struct;
                %stdvs around mean
                pccs_stdexp{iap}{im_ptr}=struct;
                pcen_stdexp{iap}{im_ptr}=struct;
                prjs_stdexp{iap}{im_ptr,ipw}=struct;
                for iar=1:2
                     lab=adj_ref_labels{iar};
                     %
                     magnifs{iap}{im_ptr,ipw}.(lab)=results_axes{id_ref,id_adj}.(lab).magnifs{im_ptr}(:,ipw);
                     %
                     switch lab
                         case 'adj'
                             z_pcc=d_adj{adj_dim};
                             z_pcc_order=sa_adj.typenames;
                         case 'ref'
                             z_pcc=d_ref{ref_dim};
                             z_pcc_order=sa_ref.typenames;
                     end
                     %
                     z_prj=results_axes{id_ref,id_adj}.(lab).projections{im_ptr}(:,:,ipw);
                     if opts_majaxes_used.plot_submean
                        z_pcc=z_pcc-repmat(mean(z_pcc,1),size(z_pcc,1),1);
                        z_prj=z_prj-repmat(mean(z_prj,1),size(z_prj,1),1);
                     end
                     for istim=1:length(z_pcc_order)
                         if length(strmatch(z_pcc_order{istim},typenames_ovlps,'exact'))~=1
                             z_pcc_order{istim}='';
                         end
                     end
                     pccs{iap}{im_ptr}.(lab)=psg_majaxes_reorder(z_pcc,stim_select,z_pcc_order);
                     coords=pccs{iap}{im_ptr}.(lab);
                     offset=mean(coords,1);
                     pcen{iap}{im_ptr}.(lab)=psg_pcaoffset(coords-repmat(offset,size(coords,1),1)); %new pc's after centering
                     prjs{iap}{im_ptr,ipw}.(lab)=psg_majaxes_reorder(z_prj,stim_select,results_axes{id_ref,id_adj}.(lab).typenames);
                     %
                     pccs_stdexp0{iap}{im_ptr}.(lab)=sqrt(mean(pccs{iap}{im_ptr}.(lab).^2,1)); %stdv around 0
                     pcen_stdexp0{iap}{im_ptr}.(lab)=sqrt(mean(pcen{iap}{im_ptr}.(lab).^2,1)); %stdv around 0
                     prjs_stdexp0{iap}{im_ptr,ipw}.(lab)=sqrt(mean(prjs{iap}{im_ptr,ipw}.(lab).^2,1)); %stdv around 0
                     %                    
                     pccs_stdexp{iap}{im_ptr}.(lab)=sqrt(var(pccs{iap}{im_ptr}.(lab),1,1)); %first 1 is to divide by N, second 1 is dimension
                     pcen_stdexp{iap}{im_ptr}.(lab)=sqrt(var(pcen{iap}{im_ptr}.(lab),1,1)); %first 1 is to divide by N, second 1 is dimension
                     prjs_stdexp{iap}{im_ptr,ipw}.(lab)=sqrt(var(prjs{iap}{im_ptr,ipw}.(lab),1,1)); %first 1 is to divide by N, second 1 is dimension
                     %
                     disp('abs value of correlations (around mean)), data pcs')
                     disp(cat(2,sprintf('%s  prj :',lab),sprintf('%8.0f',[1:adj_ref_pairs(iap,iar)])));
                     corrs=corr(pccs{iap}{im_ptr}.(lab),prjs{iap}{im_ptr,ipw}.(lab));
                     for k=1:adj_ref_pairs(iap,iar)
                         disp(cat(2,sprintf('data pc %3.0f: ',k),sprintf('%8.4f',corrs(k,:))));
                     end
                     disp('abs value of correlations (around mean)), recentered pcs')
                     disp(cat(2,sprintf('%s  prj :',lab),sprintf('%8.0f',[1:adj_ref_pairs(iap,iar)])));
                     corrs=corr(pcen{iap}{im_ptr}.(lab),prjs{iap}{im_ptr,ipw}.(lab));
                     for k=1:adj_ref_pairs(iap,iar)
                         disp(cat(2,sprintf(' rec pc %3.0f: ',k),sprintf('%8.4f',corrs(k,:))));
                     end
                end %iar
                disp(' ');
                disp(' projections onto principal components derived from data');
                for iar=1:2
                     lab=adj_ref_labels{iar};
                     v0=pccs_stdexp0{iap}{im_ptr}.(lab);
                     disp(cat(2,sprintf('%s: sqrt(var (around zero)     explained by  data pcs): ',lab),sprintf('%8.4f ',v0),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v0)/min(v0),geomean(v0)/mean(v0))));
                     v=pccs_stdexp{iap}{im_ptr}.(lab);
                     disp(cat(2,sprintf('%s: sqrt(var (around mean)     explained by  data pcs): ',lab),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                     hlid_majaxes_stats(pccs{iap}{im_ptr}.(lab));
                end
                if id_adj==id_ref
                    v=pccs_stdexp{iap}{im_ptr}.ref./pccs_stdexp{iap}{im_ptr}.adj;
                    disp(cat(2,sprintf('          sqrt(var (around mean))       ratio (ref/adj): '),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                end
                %
                disp(' ');
                disp(' projections onto principal components recomputed after centering orn and kc responses');
                for iar=1:2
                     lab=adj_ref_labels{iar};
                     v0=pcen_stdexp0{iap}{im_ptr}.(lab);
                     disp(cat(2,sprintf('%s: sqrt(var (around zero) pcs recomputed p centering): ',lab),sprintf('%8.4f ',v0),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v0)/min(v0),geomean(v0)/mean(v0))));
                     v=pcen_stdexp{iap}{im_ptr}.(lab);
                     disp(cat(2,sprintf('%s: sqrt(var (around mean) pcs recomputed p centering): ',lab),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                     hlid_majaxes_stats(pcen{iap}{im_ptr}.(lab));
                end
                if id_adj==id_ref
                    v=pcen_stdexp{iap}{im_ptr}.ref./pcen_stdexp{iap}{im_ptr}.adj;
                    disp(cat(2,sprintf('          sqrt(var (around mean))       ratio (ref/adj): '),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                end
                %
                disp(' ');
                disp(sprintf(' projections onto principal axes of transform, model component %2.0f',ipw));
                for iar=1:2
                    lab=adj_ref_labels{iar};
                    v0=prjs_stdexp0{iap}{im_ptr,ipw}.(lab);
                    disp(cat(2,sprintf('%s: sqrt(var (around zero)    explained by xform prjs): ',lab),sprintf('%8.4f ',v0),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v0)/min(v0),geomean(v0)/mean(v0))));
                    v=prjs_stdexp{iap}{im_ptr,ipw}.(lab);
                    disp(cat(2,sprintf('%s: sqrt(var (around mean)    explained by xform prjs): ',lab),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                    hlid_majaxes_stats(prjs{iap}{im_ptr,ipw}.(lab));
                end
                if id_adj==id_ref
                    v=prjs_stdexp{iap}{im_ptr,ipw}.ref./prjs_stdexp{iap}{im_ptr,ipw}.adj;
                    disp(cat(2,sprintf('          sqrt(var (around mean))       ratio (ref/adj): '),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                end
                disp(' ');
                %
                if id_adj==id_ref
                    if max(magnifs{iap}{im_ptr,ipw}.ref-magnifs{iap}{im_ptr,ipw}.adj)>tol
                        disp('warning: disagreement of magnif factors');
                    end
                end
                for iar=1:2
                    lab=adj_ref_labels{iar};
                    disp(cat(2,sprintf('%s:                       magnification factor on prjs: ',lab),sprintf('%8.4f ',magnifs{iap}{im_ptr,ipw}.(lab))));
                end
                disp(' ');
                disp('comparing axes in adj and ref space as cosines');
                disp(' projections onto principal components derived from data');
                hlid_majaxes_arc(pccs{iap}{im_ptr});
                disp(' projections onto principal components recomputed after centering orn and kc responses');
                hlid_majaxes_arc(pcen{iap}{im_ptr});
                disp(sprintf(' projections onto principal axes of transform, model component %2.0f',ipw));
                hlid_majaxes_arc(prjs{iap}{im_ptr,ipw});
                disp(' ');
             end %ipw            
             %
             disp('calculations based on model transformation from adj space')
            %compute transformation of adj to ref
            model_class=model_types_def.(model_type).class;
            im_geo=strmatch(model_type,results_geo{id_ref,id_adj}.model_types_def.model_types,'exact'); %pointer into model list in results_geo
            transform=results_geo{id_ref,id_adj}.transforms{im_geo};
            if_center=results_geo{id_ref,id_adj}.if_center;
            adj=psg_majaxes_reorder(d_adj{id_adj},stim_select,sa_adj.typenames); %select stimuli in common
            ref=psg_majaxes_reorder(d_ref{id_ref},stim_select,sa_ref.typenames);
            if if_center
                adj=adj-repmat(mean(adj,1),size(adj,1),1);
                ref=ref-repmat(mean(ref,1),size(ref,1),1);
            end
            adj_model_geo=psg_geomodels_apply(model_class,adj,transform); %modeled transformation of adj datase
            %
            %if model type is affine_offset, check that a further affine offset model will not improve
            if strcmp(model_type,'affine_offset')
                [d,adj_model_check,transform_check,opts_used]=psg_geo_affine(ref,adj_model_geo);
                disp('checking whether the affine offset model can be improved (centering ignored)');
                disp(transform_check);
                disp(transform_check.T);
            end
            adj_model_geo-repmat(mean(adj_model_geo,1),size(adj_model_geo,1),1);
            pamc{iap}{im_ptr}=psg_pcaoffset(adj_model_geo-repmat(mean(adj_model_geo,1),size(adj_model_geo,1),1)); %pcs of modeled transformation,  after centering
            pamc_stdexp0{iap}{im_ptr}=sqrt(mean(pamc{iap}{im_ptr}.^2,1)); %stdv around 0
            pamc_stdexp{iap}{im_ptr}=sqrt(var(pamc{iap}{im_ptr},1,1)); %first 1 is to divide by N, second 1 is dimension
            v0=pamc_stdexp0{iap}{im_ptr}(1:min(id_adj,id_ref));
            disp(cat(2,sprintf('modeled ref: sqrt(var (around zero)) of transformed adj: '),sprintf('%8.4f ',v0),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v0)/min(v0),geomean(v0)/mean(v0))));
            v=pamc_stdexp{iap}{im_ptr}(1:min(id_adj,id_ref));
            disp(cat(2,sprintf('modeled ref: sqrt(var (around mean)) of transformed adj: '),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
        end %im_ptr
    end %isempty
end %iap

function hlid_majaxes_stats(v)
disp(cat(2,sprintf('%57s','skewness: '),sprintf('%8.4f ',skewness(v,0,1)))); %first ar is 0 not to debias, second arg is dimension
disp(cat(2,sprintf('%57s','excess kurtosis: '),sprintf('%8.4f ',kurtosis(v,0,1)-3))); %first ar is 0 not to debias, second arg is dimension
return
end

function hlid_majaxes_arc(vecs)
disp(cat(2,'adj/ref',sprintf('%5.0f   ',[1:size(vecs.ref,2)])));
vec_a=vecs.adj./repmat(sqrt(sum(vecs.adj.^2,1)),size(vecs.adj,1),1);
vec_r=vecs.ref./repmat(sqrt(sum(vecs.ref.^2,1)),size(vecs.ref,1),1);
dots=vec_a'*vec_r;
for k=1:size(dots,1)
    disp(cat(2,sprintf('adj %1.0f',k),sprintf('%8.4f',dots(k,:))));
end
return
end


