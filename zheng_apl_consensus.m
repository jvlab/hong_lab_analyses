%zheng_apl_consensus: form consensus of TNT and control files 
%run after zheng_apl_read_align.m
%
if ~exist('opts_nonan') opts_nonan=struct(); end %for psg_remnan_coordsets
if ~exist('opts_pcon') opts_pcon=struct(); end % for procrustes_consensus
if ~exist('opts_pca') opts_pca=struct(); end % for psg_pcaoffset
%
opts_nonan=filldefault(opts_nonan,'if_log',1);
%
opts_pcon=filldefault(opts_pcon,'allow_reflection',1);
opts_pcon=filldefault(opts_pcon,'allow_offset',1);
opts_pcon=filldefault(opts_pcon,'allow_scale',0);
opts_pcon=filldefault(opts_pcon,'max_niters',1000); %nonstandard max
%
opts_pca=filldefault(opts_pca,'if_log',0);
opts_pca.nd_max=Inf;
%
%pcon_dim_max=getinp('maximum dimension for the consensus alignment dataset to be created','d',[1 max(dim_list_all)],pcon_dim_max);
%pcon_dim_max_comp=getinp('maximum dimension for component datasets to use (higher dimensions will be zero-padded)','d',[1 pcon_dim_max],pcon_dim_max);
%pcon_init_method=getinp('method to use for initialization (>0: a specific set, 0 for PCA, -1 for PCA with forced centering, -2 for PCA with forced non-centering','d',[-2 nsets],0);
pcon_dim_max=nstims_all-dim_reduce;
pcon_dim_max_comp=pcon_dim_max;
pcon_init_method=0;
%opts_pcon.allow_scale=getinp('1 to allow scaling for consensus','d',[0 1],opts_pcon.allow_scale);
disp(sprintf('pcon_dim_max=%3.0f, pcon_dim_max_comp=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',pcon_dim_max,pcon_dim_max_comp,pcon_init_method,opts_pcon.allow_scale));
opts_pcon.initialize_set=pcon_init_method;
opts_pcon.initialize_set='pca';
%
if_c2p=getinp('1 to rotate consensus into PCA space','d',[0 1]); %09Jun25
if if_c2p
    c2p_string='-pc';
else
    c2p_string='';
end
%
consensus=cell(pcon_dim_max,1);
z=cell(pcon_dim_max,1);
znew=cell(pcon_dim_max,1);
ts=cell(pcon_dim_max,1);
details=cell(pcon_dim_max,1);
opts_pcon_used=cell(pcon_dim_max,1);
%
ds_knitted=cell(1,pcon_dim_max); %reverse order of dimensions, 21Nov24
ds_components=cell(1,nsets); %partial datasets, aligned via Procrustes
%
disp('overlap matrix from stimulus matches (NaN values considered to be present')
disp(ovlp_array'*ovlp_array);
coords_isnan=zeros(nstims_all,nsets);
for iset=1:nsets
    coords_isnan(:,iset)=isnan(ds_align{iset}{1});
end
disp(sprintf('number of overlapping stimuli in component removed because coordinates are NaN'));
disp(sum(coords_isnan.*ovlp_array,1));
ovlp_array=ovlp_array.*(1-coords_isnan); %adjust overlap array to take into account NaNs (25May24)
opts_pcon.overlaps=ovlp_array;
disp('overlap matrix after excluding NaN coords in component data files')
disp(opts_pcon.overlaps'*opts_pcon.overlaps);
%
for ip=1:pcon_dim_max
    z{ip}=zeros(nstims_all,ip,nsets);
    pcon_dim_use=min(ip,pcon_dim_max_comp); %pad above pcon_dim_pad
    for iset=1:nsets
        z{ip}(:,1:pcon_dim_use,iset)=ds_align{iset}{ip}(:,[1:pcon_dim_use]); %only include data up to pcon_dim_use
        z{ip}(opts_align_used.which_common_kept(:,iset)==0,:,iset)=NaN; % pad with NaN's if no data %changed from which_common to allow for more general behavior of psg_align_coordsets when opts_align.min>1
    end
    [consensus{ip},znew{ip},ts{ip},details{ip},opts_pcon_used{ip}]=procrustes_consensus(z{ip},opts_pcon);
    disp(sprintf(' creating Procrustes consensus for dim %2.0f based on datasets up to dimension %2.0f, iterations: %4.0f, final total rms dev: %8.5f',...
        ip,pcon_dim_max_comp,length(details{ip}.rms_change),sqrt(sum(details{ip}.rms_dev(:,end).^2))));
    ds_knitted{ip}=consensus{ip};
    for iset=1:nsets
        ds_components{iset}{1,ip}=znew{ip}(:,:,iset);
    end
end
%
%implement PCA rotation if requested:  note that this is applied both to consensus{ip} and to ds_components{ip}
%
ds_knitted_orig=ds_knitted;
ds_components_orig=ds_components;
if if_c2p
    for ip=1:pcon_dim_max
        knitted_centroid=mean(ds_knitted{ip},1,'omitnan');
        [ds_knitted{ip},recon_coords,var_ex,var_tot,coord_maxdiff,opts_used_pca]=psg_pcaoffset(ds_knitted{ip},knitted_centroid,opts_pca);
%        qu=opts_used_pca.qu;
%        qs=opts_used_pca.qs;
        v=opts_used_pca.qv;
        % coords=u*s*v', and recon_coords= u*s, with v'*v=I, so recon_coords=coords*v
        for iset=1:nsets
            consensus_centroid_rep=repmat(mean(ds_components{iset}{1,ip},1,'omitnan'),nstims_all,1);
            ds_components{iset}{1,ip}=consensus_centroid_rep+(ds_components{iset}{1,ip}-consensus_centroid_rep)*v(1:ip,:);
        end
    end %ip
end


