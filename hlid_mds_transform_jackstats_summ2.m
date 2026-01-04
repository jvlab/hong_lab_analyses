%hlid_mds_transform_jackstats_summ2: summary plots for hlid_mds_transform_jackstats
% focusing on directions of principal axes, and max Euclidean dimension, jackknifed on stimuli
% (embeddings via pca, mds, cosine distances, and Pearson distances)
%
%runs on results structure from hlid_mds_transform_jackstats
%
%  See also:  HLID_MDS_TRANSFORM_JACKSTATS.
%
rng_state=rng;
if (results.if_frozen~=0) 
    rng('default');
end
colors=rand(results.nstims,3);
rng(rng_state);
%
for igp=1:results.ngps
    file_range=[1:results.nsets_gp(igp)]+sum(results.nsets_gp(1:igp-1));
    med=squeeze(min(results.max_euc_dim(file_range,:,:),[],1));
    med_jack=squeeze(min(min(results.max_euc_dim_jackstim,[],1),[],4));
    disp(sprintf(' group %1.0f: max Euclidean dimension (minimum across files)',igp))
    disp('                           data     min(jackknife by stim)')
    disp('                       sm=0  sm=1       sm=0  sm=1')
    for imeth_ptr=1:length(results.meth_use_list)
        imeth=results.meth_use_list(imeth_ptr);
        disp(sprintf('%20s  %4.0f  %4.0f       %4.0f  %4.0f',results.meth_names_short{imeth},med(imeth,:),med_jack(imeth,:)));
    end
end
%to access projections:
%results.geo_majaxes_jack_by_stim{ism,imeth,1,ijack}{idim,idim}.[ref|adj].projections{1}
%
%
for imeth_ptr=1:length(results.meth_use_list)
    imeth=results.meth_use_list(imeth_ptr);
    for isubmean_ptr=1:length(results.submean_use_list)
        isubmean=results.submean_use_list(isubmean_ptr);
        for ijack=0:results.nstims
            % if (ijack==0)
            %     mf_all=results.geo_majaxes{1+isubmean,imeth,1};
            % else
            %     mf_all=results.geo_majaxes_jack_by_stim{1+isubmean,imeth,1,ijack};
            % end
            % for k=2:results.dimlist(end)
            %     magfacs(k-1,1:2,ijack+1,1+isubmean,imeth)=mf_all{k,k}.ref.magnifs{1}(1:2)'; %take top two values
            %     magfacs(k-1,3,ijack+1,1+isubmean,imeth)=mf_all{k,k}.ref.magnifs{1}(end); %lowest value
            %     magfacs(k-1,4,ijack+1,1+isubmean,imeth)=geomean(mf_all{k,k}.ref.magnifs{1}); %lowest value
            % end
        end %ijack
    end %isubmean_ptr
end %imeth_ptr

vplot_name='angle cosines'
for ifig=1:1
    figure;
    set(gcf,'Position',[50 100 1450 800]);
    set(gcf,'Name',vplot_name);
    set(gcf,'NumberTitle','off');

    
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,vplot_name);
    axis off;
    rbase=results.geo{1+results.submean_use_list(1),results.meth_use_list(1)}{1};
    axes('Position',[0.01,0.03,0.01,0.01]); %for text
    text(0,0,sprintf('ref: %s',rbase.ref_file),'Interpreter','none');
    axis off;
    axes('Position',[0.01,0.05,0.01,0.01]); %for text
    text(0,0,sprintf('adj: %s',rbase.adj_file),'Interpreter','none');
    axis off;
end %ifig
