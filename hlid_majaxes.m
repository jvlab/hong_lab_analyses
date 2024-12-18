%hlid_majaxes: analyze the transformations of an affine geometric model to determine major axes
% for Hong Lab data
%
% needs two data files (one as reference, one to be adjusted), and the results of psg_geomodels_run that
% determined geometric models between adjusted and reference dataset
%
% 17Sep24: change defaults for plot_pairs
% 17Sep24: add options to transform principal directions to ORN space, if
%   adj data file has a field roi_names (typically from hlid_orn_merge.m), as
%   well as SVD data,  in coord_opts.aux
% 04Nov24: reorder input, so that results file is used for defaults for adj and ref files
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE, HLID_SETUP, PSG_MAJAXES, HLID_MAJAXES_COMPARE.
%
hlid_setup;
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
opts_majaxes=filldefault(opts_majaxes,'plot_pairs',[3 3;4 3;5 3;6 3;4 4]);
opts_majaxes=filldefault(opts_majaxes,'plot_order',display_orders.kcclust);
%
legacy_inputs=getinp('1 for legacy inputs','d',[0 1],0);
%
if (legacy_inputs==1)
    if getinp('1 to use hlid_orn_merge result as adj dataset','d',[0 1])
        if ~exist('fn_ref') fn_ref='./data/kc_soma_nls/hlid_odor17_coords_kc_soma_nls_consensus_scale.mat'; end
        if ~exist('fn_adj') fn_adj='./data/orn_terminals/hlid_odor17_coords_merged.mat'; end
        if ~exist('fn_geo') fn_geo='psg_geomodels_run_orn9setsMerge_kcConsensusScale_6d.mat'; end
    else
        if ~exist('fn_ref') fn_ref='./data/kc_soma_nls/hlid_odor17_coords_kc_soma_nls_consensus_scale.mat'; end
        if ~exist('fn_adj') fn_adj='./data/orn_terminals/hlid_odor17_coords_orn_terminals_consensus_scale.mat'; end
        if ~exist('fn_geo') fn_geo='psg_geomodels_run_orn7setsConsensusScale_kcConsensusScale_6d.mat'; end
        % alternate: psg_geomodels_run_orn7setsPooled_kcPooled.mat % this has piecewise affine also but is pooled, not consensus
    end
    %
    fn_ref=getinp('file name for dataset that is the reference','s',[],fn_ref);
    [d_ref,sa_ref]=psg_read_coorddata(fn_ref,[],opts_read);
    %
    fn_adj=getinp('file name for dataset that is adjusted','s',[],fn_adj);
    ds_adj=load(fn_adj);
    [d_adj,sa_adj]=psg_read_coorddata(fn_adj,[],opts_read);
    %
    fn_geo=getinp('file name containing results of psg_geomodels_run','s',[],fn_geo);
    results_geo=getfield(load(fn_geo),'results');
else
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
    ds_adj=load(fn_adj);
    [d_adj,sa_adj]=psg_read_coorddata(fn_adj,[],opts_read);
end
%
[results_axes,opts_majaxes_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts_majaxes);
%
%from hlid_coords_svd:
% SVD of responses is u*s*v', with columns of u and columns of v orthonormal
% coords_all=u*s: this is coordinates, in stimulus space
%from psg_majaxes:
% prj_dim=size(d_coords,2); %number of dimensions to project onto, less than size(A) if ref_dim<adj_dim lab='ref'
% eivecs_project=eivecs(1:prj_dim,1:prj_dim); %if id_ref<id_adj, remaining eigenvecs are units
% results_axes{id_ref,id_adj}.adj.projections{im_ptr}(:,:,ipw)=d_coords*eivecs_project;
%
% here, we want to re-express the projections in terms of coords in orn space, v*s
%
c=load(fn_adj);
if isfield(c,'roi_names') & isfield(c,'coord_opts')
    disp('determining projections in roi space')
    aux=c.coord_opts.aux;
    d_coords_roi_all=aux.v*aux.s;
    if isfield(c.coord_opts,'aux')
        aux=c.coord_opts.aux;
        nds_ref=size(results_axes,1);
        nds_adj=size(results_axes,2);
        for id_ref=1:nds_ref
            for id_adj=1:nds_adj
                if ~isempty(results_axes{id_ref,id_adj})
                    ref_dim=results_axes{id_ref,id_adj}.ref_dim;
                    adj_dim=results_axes{id_ref,id_adj}.adj_dim;
                    results_axes{id_ref,id_adj}.adj.roi_names=c.roi_names;
                    for im_ptr=1:length(results_axes{id_ref,id_adj}.model_types) %step through each model
                        d_coords_roi=d_coords_roi_all(:,1:adj_dim);
                        npw=size(results_axes{id_ref,id_adj}.adj.eivecs{im_ptr},3); %number of transformations
                        for ipw=1:npw %step through each transformation
                            eivecs=results_axes{id_ref,id_adj}.adj.eivecs{im_ptr}(:,:,ipw); %no need to adjoin extra dimensions since this is in adj space
                            results_axes{id_ref,id_adj}.adj.projections_roi{im_ptr}(:,:,ipw)=d_coords_roi*eivecs;
                        end %npw
                    end %im_ptr
                end 
            end %id_adj
        end %id_ref
        if getinp('1 to save axes information','d',[0 1])
            fn_save=getinp('file name to save axes information for comparison (e.g., hlid_majaxes_roi*.mat)','s',[]);
            opts_majaxes_used=rmfield(opts_majaxes_used,'figh');
            save(fn_save,'results_axes','fn_ref','fn_adj','fn_geo','opts_majaxes_used');
        end
    else
        disp('coord_opts.aux missing')
    end
end
