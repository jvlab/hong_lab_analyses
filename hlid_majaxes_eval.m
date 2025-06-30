%hlid_majaxes_eval: further anlaysis of transformations of an affine geometric model to determine major axes
% for Hong Lab data
%
% derived from hlid_majaxes, designed to look at transformation between ORN space and KC space, 
%   compare the major axes of the transformation with the PC's of the space
%   compare the sqrt(variance) on each PC or principal axis with the magnif factors 
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE, HLID_SETUP, PSG_MAJAXES, HLID_MAJAXES_COMPARE,
%   PSG_MAJAXES_REORDER.
%
% 30Jun25: fix bug: only consider stimuli in common to both ref and adj in computing statistics
%
hlid_setup;
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
typenames_ovlps=cell(0);
ref_ovlps=zeros(sa_ref.nstims,1);
for istim=1:sa_ref.nstims
    ref_ovlps(istim)=double(length(strmatch(sa_ref.typenames{istim},sa_adj.typenames,'exact')==1));
end
disp(sprintf('ref set has %3.0f stimuli, %3.0f in common with adj',sa_ref.nstims,sum(ref_ovlps)));
adj_ovlps=zeros(sa_adj.nstims,1);
for istim=1:sa_adj.nstims
    adj_ovlps(istim)=double(length(strmatch(sa_adj.typenames{istim},sa_ref.typenames,'exact')==1));
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
pccs=cell(n_pairs,1); %principal component coordinates
prjs=cell(n_pairs,1); %principal eigenvectors of the transformation, projected into odor space
pccs_stdexp=cell(n_pairs,1); %sqrt of variance explained
for iap=1:n_pairs
    id_adj=adj_ref_pairs(iap,1);
    id_ref=adj_ref_pairs(iap,2);
    if ~isempty(results_axes{id_ref,id_adj})
        ref_dim=results_axes{id_ref,id_adj}.ref_dim;
        adj_dim=results_axes{id_ref,id_adj}.adj_dim;
        disp(' ');
        disp(sprintf(' adj dim %2.0f ref dim %2.0f',adj_dim,ref_dim));
        magnifs{iap}=cell(0); %magnifications
        pccs{iap}=cell(0); %principal component coordinates
        prjs{iap}=cell(0); %principal eigenvectors of the transformation, projected into odor space
        pccs_stdexp{iap}=cell(0); %sqrt of variance explained
        for im_ptr=1:length(results_axes{id_ref,id_adj}.model_types) %step through each model
            npw=size(results_axes{id_ref,id_adj}.adj.eivecs{im_ptr},3); %number of transformations
            for ipw=1:npw %step through each transformation
                disp(sprintf(' model type %20s  component %1.0f',results_axes{id_ref,id_adj}.model_types{im_ptr},ipw));
                magnifs{iap}{im_ptr,ipw}=struct;
                pccs{iap}{im_ptr,ipw}=struct;
                prjs{iap}{im_ptr,ipw}=struct;
                pccs_stdexp0{iap}{im_ptr,ipw}=struct;
                prjs_stdexp0{iap}{im_ptr,ipw}=struct;
                pccs_stdexp{iap}{im_ptr,ipw}=struct;
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
                     pccs{iap}{im_ptr,ipw}.(lab)=psg_majaxes_reorder(z_pcc,stim_select,z_pcc_order);
                     prjs{iap}{im_ptr,ipw}.(lab)=psg_majaxes_reorder(z_prj,stim_select,results_axes{id_ref,id_adj}.(lab).typenames);
                     %
                     pccs_stdexp0{iap}{im_ptr,ipw}.(lab)=sqrt(mean(pccs{iap}{im_ptr,ipw}.(lab).^2,1)); %stdv around 0
                     prjs_stdexp0{iap}{im_ptr,ipw}.(lab)=sqrt(mean(prjs{iap}{im_ptr,ipw}.(lab).^2,1)); %stdv around 0
                     %                    
                     pccs_stdexp{iap}{im_ptr,ipw}.(lab)=sqrt(var(pccs{iap}{im_ptr,ipw}.(lab),1,1)); %first 1 is to divide by N, second 1 is dimension
                     prjs_stdexp{iap}{im_ptr,ipw}.(lab)=sqrt(var(prjs{iap}{im_ptr,ipw}.(lab),1,1)); %first 1 is to divide by N, second 1 is dimension
                     %
                     disp('abs value of correlations (around mean))')
                     disp(cat(2,sprintf('%s  prj :',lab),sprintf('%8.0f',[1:adj_ref_pairs(iap,iar)])));
                     corrs=corr(pccs{iap}{im_ptr,ipw}.(lab),prjs{iap}{im_ptr,ipw}.(lab));
                     for k=1:adj_ref_pairs(iap,iar)
                         disp(cat(2,sprintf('   pc %3.0f: ',k),sprintf('%8.4f',corrs(k,:))));
                     end
                end %iar
                disp(' ');
                for iar=1:2
                     lab=adj_ref_labels{iar};
                     v0=pccs_stdexp0{iap}{im_ptr,ipw}.(lab);
                     disp(cat(2,sprintf('%s: sqrt(var (around zero) explained by  pcs): ',lab),sprintf('%8.4f ',v0),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v0)/min(v0),geomean(v0)/mean(v0))));
                     v=pccs_stdexp{iap}{im_ptr,ipw}.(lab);
                     disp(cat(2,sprintf('%s: sqrt(var (around mean) explained by  pcs): ',lab),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                end
                v=pccs_stdexp{iap}{im_ptr,ipw}.ref./pccs_stdexp{iap}{im_ptr,ipw}.adj;
                disp(cat(2,sprintf('                               ratio (ref/adj): '),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                disp(' ');
                for iar=1:2
                    lab=adj_ref_labels{iar};
                    v0=prjs_stdexp0{iap}{im_ptr,ipw}.(lab);
                    disp(cat(2,sprintf('%s: sqrt(var (around zero) explained by prjs): ',lab),sprintf('%8.4f ',v0),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v0)/min(v0),geomean(v0)/mean(v0))));
                    v=prjs_stdexp{iap}{im_ptr,ipw}.(lab);
                    disp(cat(2,sprintf('%s: sqrt(var (around mean) explained by prjs): ',lab),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                end
                v=prjs_stdexp{iap}{im_ptr,ipw}.ref./prjs_stdexp{iap}{im_ptr,ipw}.adj;
                disp(cat(2,sprintf('                               ratio (ref/adj): '),sprintf('%8.4f ',v),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                disp(' ');
                if max(magnifs{iap}{im_ptr,ipw}.ref-magnifs{iap}{im_ptr,ipw}.adj)>tol
                    disp('warning: disagreement of magnif factors');
                    for iar=1:2
                        lab=adj_ref_labels{iar};
                        disp(cat(2,sprintf('%s:              magnification factor on prjs: ',lab),sprintf('%8.4f ',magnifs{iap}{im_ptr,ipw}.(lab))));
                    end
                else
                    v=magnifs{iap}{im_ptr,ipw}.(lab);
                    disp(cat(2,sprintf('                  magnification factor on prjs: '),sprintf(' max/min: %8.4f gm/am: %8.4f',max(v)/min(v),geomean(v)/mean(v))));
                end
             end %ipw
        end %im_ptr
    end %isempty
end %iap:
