function [S,afalwt_fit] = fillInNaNs(S,opts)
% JDD - fills in the NaNs 

% To use JVs code mostly as is, I need to make resps_gur using S.
numSets = length(S);
resps_gur = cell(numSets,1);
resps_tofill = cell(numSets,1);
afalwt_fit = cell(numSets,1);

for setindx = 1:numSets
    numFiles = length(S{setindx});
    [numStim,numGlom] = size(S{setindx}{1});
    resps_gur{setindx} = zeros(numStim*numGlom,numFiles);
    for fileindx = 1:numFiles
        resps_gur{setindx}(:,fileindx)=reshape(S{setindx}{fileindx}{:,:},[numStim*numGlom,1]);
    end
    
    resps_tofill{setindx} = isnan(resps_gur{setindx});
    
    afalwt_fit{setindx}=afalwt(resps_gur{setindx},1-resps_tofill{setindx},opts);
    
    resps_gur_fitted=(afalwt_fit{setindx}.x_true*afalwt_fit{setindx}.b_norm+repmat(afalwt_fit{setindx}.a,size(resps_gur{setindx},1),1)); %interpolated data
    resps_gur_filled=resps_gur{setindx};
    resps_gur_filled(resps_tofill{setindx})=resps_gur_fitted(resps_tofill{setindx});
    resps_gu_filled=reshape(resps_gur_filled,[numStim numGlom numFiles]);
    
    for fileindx = 1:numFiles
        S{setindx}{fileindx}{:,:} = resps_gu_filled(:,:,fileindx);
    end
    
    if(~opts.suppressoutput)
        fprintf('%4.0f NaN values filled in. \n',sum(resps_tofill{setindx}(:)));
    end
end


end


