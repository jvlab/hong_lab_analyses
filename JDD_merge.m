% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

opts = struct;
opts.rep_char = '+';
opts.trial_repeats = 3;
opts.set_names = {'kiTC','meTC','moSK','vaTC'};
numSets = 3;
numFiles = cell(numSets,1);

% Go to the data.
%cd ../orn_terminals_Oct25

% load the data into a structure

s = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Stot{1} = s;
numFiles{1} = length(s);

s= fileToRaw('../orn_terminals_Oct25/megamat17');
Stot{2} = s;
numFiles{2} = length(s);

s = fileToRaw('../orn_terminals_Oct25/monat');
Stot{3} = s;

% Check stimulus consistency across files within a set
[status,S,da] = checkConsist(Stot,opts);

numGlom = length(S{1,1}.Properties.VariableNames);

missing_glomeruli = cell(numSets,1);

for setindx = 1:length(da)
    
    missing_glomeruli{setindx} = zeros(numGlom,length(da{setindx}));
    for fileindx = 1:length(da{setindx})
        nanCols = find(all(isnan(S{setindx,fileindx}{:,:})));
        missing_glomeruli{setindx}(nanCols,fileindx)=1;        
    end
end

for setindx = 1:numSets
    numFiles = length(da{setindx});
    fprintf('Set %s \n',opts.set_names{setindx});
    for fileCount = 0:numFiles    
        fprintf('number of glomeruli present in at least %1.0f datasets: %1.0f (missing in %1.0f or less) \n',fileCount,sum(sum(missing_glomeruli{setindx},2)<=numFiles-fileCount),numFiles-fileCount);
    end
end
    






