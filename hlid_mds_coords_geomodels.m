%hlid_mds_coords_geomodels:  examine geometric models for transformations
%between hlid datasets in which spaces have been constructed via alternative multidimensional scaling methods
% (using hlid_rastim_mds_coords_demo)
%
% variants include:
%   Euclidean distances between z-scored values, with dimension reduction
%     via svd of z-scored values, or mds of distances
%   angular distance (arccos) between dot products of normalized z-scored values 
%   arc length distance between normalized z-scored values
%   angular distance (arccos) between dot products of normalized z-scored values after subtraction of mean across ROIs (Pearson)
%   arc length distance between normalized z-scored values ater subtraction of mean across ROIs (Pearson)
%
% uses the workspaces saved by hlid_rastim_mds_coords_demo)
%
% 09Dec25: option to shnorten  r{isub,imeth}.results by removing opts_model_shuff_used_nestdim
%
%  See also:  HLID_LOCALOPTS, HLID_RASTIM_MDS_COORDS_DEMO, PSG_GEOMODELS_RUN.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
%dataset 1 will be reference, dataset 2 will be adjusted to fit.
%
if_debug=getinp('1 for debugging','d',[0 1]);
if (if_debug)
    d_max=5;
    max_iters=10;
    nshuffs_def=10;
    ref_dim_list_def=[2:3];
else
    d_max=10;
    max_iters=1000;
    nshuffs_def=100;
    ref_dim_list_def=[1:7];
end
shorten_field='opts_model_shuff_used_nestdim';
if_shorten=getinp(sprintf('shorten results structure by removing ''%s''',shorten_field),'d',[0 1],1);
ra_strings={'ref','adj'};
sm_strings={'nosub mean','sub mean'};
nra=length(ra_strings);
if ~exist('ref_adj_filenames')
    ref_adj_filenames=cell(nra,1);
    ref_adj_filenames{1}='.\plots\hlid_rastim_mds_kc-tnt3c_07Nov25.mat';
    ref_adj_filenames{2}='.\plots\hlid_rastim_mds_kc-tntlabel_07Nov25.mat';
end
ndims_avail=Inf;
dra=cell(nra,1);
for ira=1:nra
    ref_adj_filenames{ira}=getinp(sprintf('mat-file with alternative embeddings for %s',ra_strings{ira}),'s',[],strrep(ref_adj_filenames{ira},'\','/'));
    dra{ira}=load(ref_adj_filenames{ira});
    %assume thesea are consistent across ref and adj
    nstims=dra{ira}.r.nstims;
    nmeths=length(dra{ira}.r.meths);
    nsubs=size(dra{ira}.r.data_knit,2);
    if_submean=nsubs-1;
    disp(sprintf('loaded with %2.0f methods, %2.0f stimuli, %2.0f original files',...
        nmeths,nstims,length(dra{ira}.r.filenames_short)));
    dra{ira}.r=rmfield(dra{ira}.r,'knit_stats');
    disp('knit_stats field removed');
    disp('setup for embedding from raw data:');
    disp(dra{ira}.knit_stats_setup);
    for imeth=1:nmeths
        for isub=1:nsubs
            ndims_avail=min(ndims_avail,length(dra{ira}.r.data_knit{imeth,isub}.ds{1}));
        end
    end
end
d_max=getinp('maximum dimension to process','d',[2 ndims_avail],min(d_max,ndims_avail));
if ~exist('neigs_show_max') neigs_show_max=10; end
%
if ~exist('opts_geofit')
    opts_geofit=struct;
end
opts_geofit.model_types_def=psg_geomodels_define(1);
opts_geofit.nshuffs=getinp('number of shuffles for geometric model fitting','d',[0 1000],nshuffs_def);
opts_geofit.ref_dim_list=getinp('list of dimensions to analyze','d',[1 d_max],ref_dim_list_def);
opts_geofit.if_log=getinp('to show logging for geometric model fitting','d',[0 1],1);
opts_geofit.if_summary=getinp('1 to show summary for geometric model fitting','d',[0 1],1);
opts_geofit.if_nestbydim=-1; %nonstandard, to allow for nesting
%
%now align to a consensus across embeddings, for ref and adj
%
dnames={'ds','sas','sets'};
aux_knit=struct;
aux_knit.opts_knit.if_pca=1; %rotate into pca space
aux_knit.opts_knit.dim_max_in=d_max;
aux_knit.opts_knit.max_iters=max_iters;
aux_knit.opts_knit.if_log=0; %nondefault
%set up plot
aux_disp=struct;
aux_disp.opts_disp.if_legend=1;
for k=1:nmeths
    aux_disp.opts_disp.set_labels{k}=dra{1}.r.meths{k}.name_short;
end
aux_disp.opts_disp.set_colors={'k','k','r','r','r','g','g','g'};
aux_disp.opts_disp.set_markers={'o','.','x','<','s','x','<','s'};
aux_disp.opts_disp.connect_sets_method='list';
aux_disp.opts_disp.connect_sets_list=[1 2;2 5;2 8;3 4;4 5;5 3;6 7;7 8;8 6;5 8]; %connect among families, connect chord with chord
aux_disp.opts_disp.connect_sets_linestyles={':'};
aux_disp.opts_disp.dim_select=4;
aux_disp.opts_disp.fig_position=[100 100 1400 900];
aux.disp.opts_disp.callout_amount=0.5;
%
components=cell(nra,nsubs);
if_ok=1;
for ira=1:nra
    for isub=1:nsubs
        %create a stack of datasets, one for each embedding
        allembeds=struct;
        disp(sprintf('assembling and checking %2.0f embedding methods for %s, %s',nmeths,ra_strings{ira},sm_strings{isub}));
        for id=1:length(dnames)
            dname=dnames{id};
            allembeds.(dname)=cell(nmeths,1);
            for imeth=1:nmeths
                allembeds.(dname){imeth}=dra{ira}.r.data_knit{imeth,isub}.(dname){1};
            end
        end
        %adjust metadata fields 
        for imeth=1:nmeths
            allembeds.sets{imeth}.type='data'; % output file did not have had a type field
            allembeds.sets{imeth}.label=sprintf('consensus for method %s',dra{ira}.r.meths{imeth}.name_full);
            allembeds.sets{imeth}.dim_list=[1:length(dra{ira}.r.data_knit{imeth,isub}.ds{1})];
            allembeds.sets{imeth}.label_long=allembeds.sets{imeth}.label;
            allembeds.sets{imeth}.pipeline=[];
            for istim=1:nstims
                typename=allembeds.sas{imeth}.typenames{istim};
                firstblank=min(find(typename==' '));
                if ~isempty(firstblank)
                    typename=typename(1:firstblank-1);
                end
                allembeds.sas{imeth}.typenames{istim}=typename;
            end
        end
        [check,opts_used]=rs_check_coordsets(allembeds,setfield(struct,'if_warn',1));
        if ~isempty(check.warnings)
            if_ok=0;
        end
        if (if_ok)
            [knitted,aux_knit_out]=rs_knit_coordsets(allembeds,aux_knit);
            components{ira,isub}=aux_knit_out.components;
            disp('knitted without re-alignment, and components extracted.')
            clear aux_knit_out;
            %quick plot
            aux_disp.opts_disp.fig_name=sprintf('%s, %s',ra_strings{ira},sm_strings{isub});
            aux_disp_out=rs_disp_coordsets(components{ira,isub},aux_disp);
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,aux_disp.opts_disp.fig_name,'FontSize',8);
            axis off;
        end
    end %isub
end %ira
%now, look at geometric model for each embedding, and each mean-subtraction
r=cell(nsubs,nmeths);
for isub=1:nsubs
    for imeth=1:nmeths
        method_desc=sprintf('%s, %s',dra{1}.r.meths{imeth}.name_full,sm_strings{isub});
        disp('**************');
        disp(sprintf('geometric modeling for %s',method_desc));
        r{isub,imeth}=struct;
        r{isub,imeth}.ref_adj_filenames=ref_adj_filenames;
        r{isub,imeth}.method_desc=method_desc;
        d_ref=components{1,isub}.ds{imeth};
        d_adj=components{2,isub}.ds{imeth};
        [results_keep,opts_geofit_used]=psg_geomodels_fit(d_adj,d_ref,opts_geofit);
        if if_shorten
            if isfield(results_keep,shorten_field)
                results_keep=rmfield(results_keep,shorten_field);
            end
        end
        r{isub,imeth}.results=results_keep;
        if ~isempty(opts_geofit_used.warnings)
            disp(opts_geofit_used.warnings);
        end
    end
end
