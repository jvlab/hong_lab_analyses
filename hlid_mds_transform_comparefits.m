%hlid_mds_transform_comparefits: demonstration of comparison of fits
%of various embeddings across model dimensions
%and shows how to retrieve quantitative data from analysis files
%
%  See also: HLID_MDS_TRANSFORM_STATS, HLID_MDS_COORDS_GEOMODELS.
load plots\hlid_mds_transform_stats_01Jan26.mat;
whos
disp('results');
disp(results);
disp('results_geo');
disp(results.geo);
disp('resutls.geo_majaxes');
disp(results.geo_majaxes);
%
dim_max=length(results.geo{1,1});
nembeds=size(results.geo,2);
nsm=size(results.geo,1);
embed_names=results.meth_names_short;
model_name=results.model_types_use.model_types{1}; %only one model
dvals=zeros(nembeds,dim_max,nsm);
for ism=1:nsm
    switch ism
        case 1
            disp('mean not subtracted');
        case 2
            disp('mean subtracted');
    end
    for iembed=1:nembeds
        rg=results.geo{ism,iembed};
        for idim=1:dim_max
            dvals(iembed,idim,ism)=rg{idim,idim}.d; %only one model
        end
        disp(sprintf('%25s %s',embed_names{iembed},sprintf('%7.4f',dvals(iembed,:,ism))));
    end
end
