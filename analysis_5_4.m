% Analysis done may 5th, 2026

Xraw = load('april20_merged_data.mat');
Xinter = Xraw.merged_data{1};

repeat_stim = findRepeatStimuli(Xinter);

LHS = Xinter(repeat_stim(1:2:length(repeat_stim)),:);
RHS = Xinter(repeat_stim(2:2:length(repeat_stim)),:);

