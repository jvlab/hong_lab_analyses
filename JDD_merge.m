% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

% dependencies : fileToRaw.m, checkConsist.m

opts = struct;
opts.rep_char = '+';
opts.trial_repeats = 3;
opts.set_names = {'kiTC','meTC','vaTC'};
opts.kdf = {'max_peak','mean_peak'}; % Known data fields.
opts.suppressoutput = true;
opts.interactive = false;
opts.restore_size = true;
opts.submean = false;
opts.hist_quantiles = [0.05 .25 .5 .75 .96];
opts.hist_bins = 50;

% I am assuming that the files that are part of a set are in a separate
% directory. There is not name checking, if a data file is in the folder
% the code will attempt to load it into the set.
Sraw{1} = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Sraw{2} = fileToRaw('../orn_terminals_Oct25/megamat17');
%Sraw{3} = fileToRaw('../orn_terminals_Oct25/monat');
Sraw{3} = fileToRaw('../orn_terminals_Oct25/validation2');

% Check stimulus consistency across files within a set,
% and glomerus consistency across sets.
[~,Sall] = checkConsist(Sraw,opts);

% First, I want to check that each stimulus has a non-NaN value somewhere
% within a set. If an all NaN is found, that stimulus is removed.
% The presence of the glomeruli across sets is determined. the glomeruli
% used appear in a minimum number of files. 
Strimmed = lookForSetWideHoles(Sall,opts);

% Fill in the remaining holes.
% The calls the afalwt interpolator.
[Sfilled,afalwt_fit] = fillInNaNs(Strimmed,opts);

% Generate the first set of plots (raw - trimmed - filled)
makePlots_1(Sall,Strimmed,Sfilled);

%
% also create the resps_set cell array
[resps_set,resp_range] = calcResp(Sfilled,afalwt_fit,opts);

% Generate the quantile plots
makePlots_quantile(resps_set,resp_range,opts);

% Since the odorants within a set are unique and each set is uniquely
% tagged, we don't need to do all of the uniqueness checks on the stimuli.
% I do need to assemble the odorants into a single usnique structure.


