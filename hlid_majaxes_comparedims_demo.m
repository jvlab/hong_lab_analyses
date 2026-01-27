%hlid_majaxes_comparedims_demo: demonstration of comparison of major axis analysis across model dimensions
%also shows how to retrieve quantitative data from analysis files
%
load hlid_majaxes_TNT_label_TNT3c_ConsensusNoScale_axes;
load psg_geomodels_run_TNT_label_TNT3c_ConsensusNoScale;
whos
disp('results');
disp(results{1,1});
disp('results_axes');
disp(results_axes{1,1});
ar={'adj','ref'};
nar=length(ar);
dim_max=length(results);
mname='affine_offset';
rm=strmatch(mname,results{1,1}.model_types_def.model_types,'exact'); %pointer to model in results
ra=strmatch(mname,results_axes{1,1}.model_types,'exact'); %pointer to model in axes
%
disp('stimulus order')
disp(results_axes{1,1}.adj.typenames);
nstims=length(results_axes{1,1}.adj.typenames);
%
projs_norm=cell(dim_max,nar); %normalized projections
for id=1:dim_max
    adj_dim=id;
    ref_dim=id;
    disp(' ');
    disp(sprintf('dimension %1.0f: transformation from adj to ref, coords are pcs in each space',id));
    T=results{adj_dim,ref_dim}.transforms{rm}.T;
    disp(T);
    for iar=1:nar
        disp(sprintf('dot-products of projections with themselves in %s',ar{iar}));
        projs=results_axes{adj_dim,ref_dim}.(ar{iar}).projections{ra};
        projs_norm{id,iar}=projs./repmat(sqrt(sum(projs.^2,1)),nstims,1);
        disp(projs_norm{id,iar}'*projs_norm{id,iar});
        %
        if (id>=2)
            disp(sprintf('dot-products of projections in %s with lower-dim %s space',ar{iar},ar{iar}))
            disp(projs_norm{id,iar}'*projs_norm{id-1,iar});
        end
    end
    disp('dot-products of projections in adj space (rows) with projections in ref space (cols)');
    disp(projs_norm{id,1}'*projs_norm{id,2});
end
