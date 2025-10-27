% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

opts = struct;
opts.rep_char = '+';

numSets = 2;
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

% Check stimulus consistency across files within a set
[status,S] = checkConsist(Stot);


    






