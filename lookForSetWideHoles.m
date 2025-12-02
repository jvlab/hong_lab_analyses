function S = lookForSetWideHoles(S,opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numSets = length(S);

% Remove all stimuli that are NaN for every set.
% I think that I want this after the next part, since removal of glom can
% result in something resolved here being unresolved.



% Glomerus stuff
for setindx = 1:numSets
    numFiles = length(S{setindx});
    nGlomeruli = length(S{setindx}{1}.Properties.VariableNames);
    missing_glomeruli = zeros(nGlomeruli,numFiles);
    
    for fileindx = 1:numFiles
        nanCols = all(isnan(S{setindx}{fileindx}{:,:}),1);
        missing_glomeruli(nanCols,fileindx) = 1;
    end
    if(~opts.suppressoutput)
        for ipresent = 0:numFiles
            fprintf('number of glomeruli present in at least %1.0f datasets: %1.0f (missing in %1.0f or less) \n',...
                ipresent,sum(sum(missing_glomeruli,2)<=numFiles-ipresent),numFiles-ipresent);
        end
    end
    % Here is where the number of files can be chosen. It works, I just use
    % a value of 1.
    if(opts.interactive)
        min_present = getinp('minimum number of preps that a glomerulus must be present in','d',[0 numFiles],1);
    else
        min_present = 1;
    end
    
    for fileindx = 1:numFiles
        S{setindx}{fileindx}=S{setindx}{fileindx}(:,sum(missing_glomeruli,2)<=numFiles-min_present);
    end
    
end

% These should be two separate functions - so the order can be switched.


for setindx = 1:numSets
    numFiles = length(S{setindx});
    [numStim,numGlom] = size(S{setindx}{1});
    bigMat = zeros(numStim*numGlom,numFiles);
    for fileindx=1:numFiles
        bigMat(:,fileindx) = reshape(S{setindx}{fileindx}{:,:},[numStim*numGlom, 1]);
    end
    
    to_remove = all(isnan(bigMat),2);
    
    to_remove_reshape = reshape(to_remove,[numStim, numGlom]);
    
    to_last = any(to_remove_reshape,2);
    
    for fileindx = 1:numFiles
        S{setindx}{fileindx}(to_last,:) = [];
    end
end


%{
for setindx = 1:numSets
    numFiles = length(S{setindx});
    [numStim,numGlom] = size(S{setindx}{1});
    totalVec = ones(numStim,1);
    for fileindx = 1:numFiles
        tmp = all(isnan(S{setindx}{fileindx}{:,:}),2);
        totalVec = totalVec.*tmp;
        
    end
    if(sum(totalVec)>0)        
        for fileindx = 1:numFiles
            S{setindx}{fileindx}(logical(totalVec),:) = [];
        end
    end
end
%}

end