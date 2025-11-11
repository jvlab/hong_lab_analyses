function [resps_set,resp_range] = calcResp(S,afalwt_fit,opts)
% JDD 11-3
% Calculates the responses and bounds

numSets = length(S);

resps_set = cell(numSets,1);
resp_range = [Inf -Inf];
for setindx = 1:numSets
    resps_set{setindx} = S{setindx}{1};
    numFiles = length(S{setindx});
    [numStim,numGlom] = size(S{setindx}{1});
    resps_set{setindx}{:,:}=reshape(afalwt_fit{setindx}.x_true,[numStim numGlom]); %use regression slope as response measure

    if opts.restore_size %restorr size if needed
        resps_set{setindx}{:,:}=resps_set{setindx}{:,:}*geomean(afalwt_fit{setindx}.b_norm);
    end
    
    resp_range(1) = min(resp_range(1),min(min(resps_set{setindx}{:,:})));
    resp_range(2) = max(resp_range(2),max(max(resps_set{setindx}{:,:})));
%
%condition the data and stimulus names
%
    if (opts.submean)
        resps_set{setindx}{:,:} = resps_set{setindx}{:,:}-repmat(mean(resps_set{setindx}{:,:},1),numStim,1);
    end 
    
    if(~opts.suppressoutput)
        fprintf('set %2.0f has %3.0f files used (of %3.0f), %3.0f stimuli, and %3.0f glomeruli used \n',...
            setindx,numFiles,numFiles,numStim,numGlom);
    end
end

end

