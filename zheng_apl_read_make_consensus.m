%zheng_apl_read_make_consensus: read, and align TNT and control files and make consensus
% no option for pca rotation, scaling, or normalization
%
%   See also:  RS_GEOFIT.
%
%read a set to use as the adjusted dataset
data_uses.adj='TNT_label';
data_uses.ref='TNT3c';
%
data_use=data_uses.adj;
zheng_apl_read_align;
d.adj.ds=ds_align;
d.adj.sas=sas_align;
d.adj.sets=sets_align;
clear sets_align sas_align ds_align 
%read a set to use as the reference
data_use=data_uses.ref;
zheng_apl_read_align;
d.ref.ds=ds_align;
d.ref.sas=sas_align;
d.ref.sets=sets_align;
clear sets_align sas_align ds_align 
%
%calculate consensus for adj and ref sets
%
if ~exist('opts_pcon') opts_pcon=struct(); end % for procrustes_consensus
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0);
opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
if ~exist('opts_pca') opts_pca=struct(); end % for psg_pcaoffset
opts_pca=filldefault(opts_pca,'if_log',0);
opts_pca.nd_max=Inf;
%
%pcon_dim_max=getinp('maximum dimension for the consensus alignment dataset to be created','d',[1 max(dim_list_all)],pcon_dim_max);
%pcon_dim_max_comp=getinp('maximum dimension for component datasets to use (higher dimensions will be zero-padded)','d',[1 pcon_dim_max],pcon_dim_max);
%pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering','d',[-2 nsets],0);
pcon_dim_max=nstims_all-dim_reduce;
pcon_dim_max_comp=pcon_dim_max;
pcon_init_method=0;
disp(sprintf('pcon_dim_max=%3.0f, pcon_dim_max_comp=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',pcon_dim_max,pcon_dim_max_comp,pcon_init_method,opts_pcon.allow_scale));
opts_pcon.initialize_set=pcon_init_method;
opts_pcon.initialize_set='pca';
%
consensus=cell(1,pcon_dim_max);
z=cell(pcon_dim_max,1);
cons=struct;
%
for iar=1:2
    switch iar
        case 1
            iar_string='adj';
        case 2
            iar_string='ref';
    end
    nsets=length(d.(iar_string).ds);
    disp(sprintf('forming consensus for %s, %2.0f datasets',iar_string,nsets));
    for ip=1:pcon_dim_max
        z{ip}=zeros(nstims_all,ip,nsets);
        pcon_dim_use=min(ip,pcon_dim_max_comp); %pad above pcon_dim_pad
        for iset=1:nsets
            z{ip}(:,1:pcon_dim_use,iset)=d.(iar_string).ds{iset}{ip}(:,[1:pcon_dim_use]); %only include data up to pcon_dim_use
        end
        [consensus{1,ip},znew,ts,details,opts_pcon_used{ip}]=procrustes_consensus(z{ip},opts_pcon);
        disp(sprintf(' creating Procrustes consensus for dim %2.0f based on datasets up to dimension %2.0f, iterations: %4.0f, final total rms dev: %8.5f',...
            ip,pcon_dim_max_comp,length(details.rms_change),sqrt(sum(details.rms_dev(:,end).^2))));
        clear znew ts details
    end
    cons.(iar_string)=consensus;
end %iar
%%%%%%%%%%%%%%%%split out from here down to allow for bootstrap analyses,
%%%%%%%%%%%%%%%%projections of principal axes
%%%%ask about saving mat file and figures
%
%now construct geometric models
%
if ~exist('model_dim_max') model_dim_max=6; end
if (if_debug)
    if ~exist('nshuffs') nshuffs=100; end
else
    if ~exist('nshuffs') nshuffs=1000; end
end
%
model_list={'procrustes_noscale_offset'  'procrustes_scale_offset'  'affine_offset'};
data_in=struct;
data_in.ds{1}=cons.adj;
data_in.sas{1}=d.adj.sas{1};
data_in.sets{1}=d.adj.sets{1};
%
data_out=struct;
data_out.ds{1}=cons.ref;
data_out.sas{1}=d.ref.sas{1};
data_out.sets{1}=d.ref.sets{1};
%
aux=struct;
%
opts_geof=struct;
opts_geof.nshuffs=nshuffs;
opts_geof.model_list=model_list;
opts_geof.dimpairs_list=repmat([1:model_dim_max]',1,2);
opts_geof.dim_max_in=model_dim_max;
opts_geof.dim_max_out=model_dim_max;
opts_geof.dimpairs_method='equal'; %only diagonal values
opts_geof.if_fit_summary=0;
opts_geof.if_log=0;
%
[gfs,xs,aux_out]=rs_geofit(data_in,data_out,setfield(aux,'opts_geof',opts_geof));
%
opts_dgeo=struct;
opts_dgeo.colors_models={'k','c','m'};
if ~exist('if_showquant') if_showquant=0; end %set to show the 0.05 quantile
opts_dgeo.if_showquant=if_showquant;
%
aux_out_disp=rs_disp_geofit(gfs{1}.gf,setfield(aux,'opts_dgeo',opts_dgeo));
%
disp(sprintf('Model comparison via Procrustes d, nshuffs=%4.0f, if_sm=%1.0f',nshuffs,if_sm));
disp(data_uses);
disp('    dim     d(affine)  p(affine v unif sc)   d(unif sc) p(unif sc v no sc)    d(no sc)');
disp('                        orig den  shuff den              orig den  shuff den');
ptr_affine=strmatch('affine_offset',model_list,'exact');
ptr_uniscale=strmatch('procrustes_scale_offset',model_list,'exact');
ptr_noscale=strmatch('procrustes_noscale_offset',model_list,'exact');
for k=1:model_dim_max
    d_shuff_affine_vs_uniscale=squeeze(gfs{1}.gf{k,k}.d_shuff(ptr_affine,:,ptr_uniscale,:));
    d_shuff_uniscale_vs_noscale=squeeze(gfs{1}.gf{k,k}.d_shuff(ptr_uniscale,:,ptr_noscale,:)); 
    d_affine=gfs{1}.gf{k,k}.d(ptr_affine);
    d_uniscale=gfs{1}.gf{k,k}.d(ptr_uniscale);
    d_noscale=gfs{1}.gf{k,k}.d(ptr_noscale);
    vars=[k d_affine sum(double(d_affine>d_shuff_affine_vs_uniscale))/nshuffs d_uniscale sum(double(d_uniscale>d_shuff_uniscale_vs_noscale))/nshuffs d_noscale];
    disp(sprintf('%8.4f   ',vars))
end
