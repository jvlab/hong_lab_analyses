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
% 
Stot{1} = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Stot{2} = fileToRaw('../orn_terminals_Oct25/megamat17');
Stot{3} = fileToRaw('../orn_terminals_Oct25/monat');
Stot{4} = fileToRaw('../orn_terminals_Oct25/validation2');

% Check stimulus consistency across files within a set,
% and glomerus consistency across sets.
[~,S] = checkConsist(Stot,opts);

% First, I want to check that each stimulus has a non-NaN value somewhere
% within a set. If found, that stimulus is removed.

S = lookForSetWideHoles(S,opts);









