function S = fillInNaNs(S,opts)
% JDD - fills in the NaNs 

% To use JVs code mostly as is, I need to make resps_gur using S.
numSets = length(S);
resps_raw = cell(numSets,1);

for setindx = 1:numSets
    numFiles = length(S{setindx});
    [numStim,numGlom] = size(S{setindx}{1});
    resps_raw{setindx} = zeros(numStim,numGlom,numFiles);
    for fileindx = 1:numSets
        resps_raw{setindx}(:,:,fileindx)=S{setindx}{fileindx}{:,:};
    end
end



%[afalwt_fit,afalwt_b_change,afalwt_optsused]=afalwt(resps_gur,1-resps_tofill,afalwt_opts);
%resps_gur_fitted=(afalwt_fit.x_true*afalwt_fit.b_norm+repmat(afalwt_fit.a,size(resps_gur,1),1)); %interpolated data
%resps_gur_filled=resps_gur;
%resps_gur_filled(resps_tofill)=resps_gur_fitted(resps_tofill);
%resps_gu_filled=reshape(resps_gur_filled,[nstims(iset) nglomeruli_use(iset) nfiles_use(iset)]);


end

