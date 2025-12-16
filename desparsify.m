function [Sn,fileList] = desparsify(S)
% JDD 12/2/25
% Written specifically to make monat data set useable. 

% S is a set of tables corresponding to the files: S{1} is the first file..

numFiles = length(S);

% Find the the number of files per odor

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

% This sets the number of synthetic files to the maximum.
% Smaller numbers can be chosen. 
% For example:
% numSynFiles = min(filecount);
% Would provide the number of files needed so that all odors are
% represented in each file. 
%fileList
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


end

