function [pattern,pcount] = desparsify(S)
% JDD 12/2/25
% Written specifically to make monat data set useable. 

% S is a set of tables corresponding to the files: S{1} is the first file..

tableList = {};

numUniqueSets = 0;
pattern = {};
pcount = {};
numFiles = length(S);
fileList = {};


for fileindx=1:numFiles
    pattern_it = all(isnan(S{fileindx}{:,:}),2);
    
    patseen = false;
    for patternIndx = 1:length(pattern)
        
        if(pattern_it==pattern{patternIndx})
            % The pattern has been seen
            pcount{patternIndx} = pcount{patternIndx}+1;
            patseen = true;
            continue;
        end
    end
    if(~patseen)
        pattern = [pattern, pattern_it];
        pcount = [pcount, 1];
    end
    
end




end

