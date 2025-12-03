function Sn = desparsify(S)
% JDD 12/2/25
% Written specifically to make monat data set useable. 

% S is a set of tables corresponding to the files: S{1} is the first file..

%pattern = {};
%pcount = {};
numFiles = length(S);


%{
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
%}
% Look at the number of files per odor
[numStim,numGlom] = size(S{1});
filecount = zeros(numStim,1);
fileList = cell(numStim,1);

for fileindx = 1:numFiles
    pattern_it = ~all(isnan(S{fileindx}{:,:}),2);
    for stimindx=1:numStim
        if(pattern_it(stimindx))
            fileList{stimindx} = [fileList{stimindx} fileindx];
            filecount(stimindx) = filecount(stimindx) + 1;
        end
    end
end

numSynFiles = max(filecount);

Sn = cell(numSynFiles,1);

for fileindx = 1:numSynFiles
    tmpArray = NaN(numStim,numGlom);
    Sn{fileindx} = array2table(tmpArray);
    Sn{fileindx}.Properties.VariableNames = S{1}.Properties.VariableNames;
    Sn{fileindx}.Properties.RowNames = S{1}.Properties.RowNames;
    
    for stimindx=1:numStim
        if(fileindx<=filecount(stimindx))
            Sn{fileindx}(stimindx,:) = S{fileList{stimindx}(fileindx)}(stimindx,:);
        end
    end
end



    


% Every odor is accounted for at least six times.


end

