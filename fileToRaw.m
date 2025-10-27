function S = fileToRaw(directoryName)
% JDD
% Loads the files in the directory.


% Set up the system command to identify and load the data files.
% It might be obvious I struggled with the matlab string vs char.

% List the files in the directory.
command = string('ls ');

command = command+string(directoryName)+string('/*.mat');

[status,fileList] = system(char(command));

if(status)
    error('listing of directory contents failed.');
end

fileList = string(fileList);

fileArray = split(fileList);

fileArray(end) = [];

numFiles = length(fileArray);

S = cell(numFiles,1);

for file = 1:numFiles
    loadString = fileArray(file);
    load(char(loadString));
    S{file}.description = description;
    S{file}.meta = meta;
    S{file}.movie_timeseries = movie_timeseries;
    S{file}.response_amplitude_stim = response_amplitude_stim;
    S{file}.response_amplitude_trials = response_amplitude_trials;
    S{file}.rois = rois;
    S{file}.trial_info = trial_info;
    S{file}.trial_timeseries = trial_timeseries;
end

end

