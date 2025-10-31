% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

% dependencies : fileToRaw.m, checkConsist.m

opts = struct;
opts.rep_char = '+';
opts.trial_repeats = 3;
opts.set_names = {'kiTC','meTC','moSK','vaTC'};
opts.kdf = {'max_peak','mean_peak'}; % Known data fields.
opts.suppressoutput = true;
opts.interactive = false;
% load the data into a structure
% Set by set.
% I am assuming that the files that are part of a set are in a separate
% directory. There is not name checking, if a data file is in the folder
% the code will attempt to load it into the set.
Sraw{1} = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Sraw{2} = fileToRaw('../orn_terminals_Oct25/megamat17');
Sraw{3} = fileToRaw('../orn_terminals_Oct25/monat');
Sraw{4} = fileToRaw('../orn_terminals_Oct25/validation2');

% Check stimulus consistency across files within a set,
% and glomerus consistency across sets.
[~,S] = checkConsist(Sraw,opts);

% First, I want to check that each stimulus has a non-NaN value somewhere
% within a set. If found, that stimulus is removed.
% The presence of the glomeruli across sets is determined. the glomeruli
% used appear in a minimum number of files. 
S = lookForSetWideHoles(S,opts);

% Fill in the remaining holes.
% The calls the afalwt commands
S = fillInNaNs(S,opts);









