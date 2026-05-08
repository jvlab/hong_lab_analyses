% Normalize a merged set
load('april20_merged_data.mat');

X = merged_data{1}; % I don't like the union because of the zeros.

% Divide this into the constituent sets.
[numOdors,numGlom] = size(X);
suffix = cell(numOdors,1);
for indx = 1:numOdors
    suffix{indx} = X.Properties.RowNames{indx}(end-3:end);
end
setID = unique(suffix);
setX = cell(length(setID),1);
for indx = 1:length(setID)
    S = contains(suffix,setID(indx));
    setX{indx} = X(S,:);
end

% Set by set we want to set the variance of the output to be the mean
% variance of the input.
norm_setX = cell(length(setID),1);
glom_variance = zeros(length(setID),numGlom);

for indx=1:length(setX)
    Xarray = setX{indx}{:,:};    
    Xarray_centered = Xarray - repmat(mean(Xarray),[size(Xarray,1) 1]);
    glom_variance(indx,:) = diag(Xarray_centered'*Xarray_centered)'/size(Xarray,1);
end
glom_mov = mean(glom_variance);

for indx = 1:length(setX)
    Xarray = setX{indx}{:,:};    
    Xarray_centered = Xarray - repmat(mean(Xarray),[size(Xarray,1) 1]);
    Xarray_normalized = sqrt(glom_mov).*(sqrt(size(Xarray,1))*Xarray_centered./repmat(sqrt(diag(Xarray_centered'*Xarray_centered))',[size(Xarray_centered,1) 1]));    
    norm_setX{indx} = setX{indx};
    norm_setX{indx}{:,:} = Xarray_normalized;
end

merged_normalize = table;
merged_normalize = norm_setX{1};
for indx = 2:length(setID)
    merged_normalize = [merged_normalize ; norm_setX{indx}];
end


% Now the diagnostic odors
load('april21_merged_data.mat');
Xd = merged_data{1}; % I don't like the union because of the zeros.

Xd = Xd(findRepeatStimuli(Xd),:);
% Divide this into the constituent sets.
[numOdors,numGlom] = size(Xd);
suffix = cell(numOdors,1);
for indx = 1:numOdors
    suffix{indx} = Xd.Properties.RowNames{indx}(end-3:end);
end
setID = unique(suffix);
setXd = cell(length(setID),1);
for indx = 1:length(setID)
    S = contains(suffix,setID(indx));
    setXd{indx} = Xd(S,:);
end

% Set by set we want to set the variance of the output to be the mean
% variance of the input.
norm_setXd = cell(length(setID),1);
glom_varianced = zeros(length(setID),numGlom);

for indx=1:length(setX)
    Xarrayd = setXd{indx}{:,:};    
    Xarray_centeredd = Xarrayd - repmat(mean(Xarrayd),[size(Xarrayd,1) 1]);
    glom_varianced(indx,:) = diag(Xarray_centeredd'*Xarray_centeredd)'/size(Xarrayd,1);
end
glom_movd = mean(glom_varianced);

for indx = 1:length(setX)
    Xarrayd = setXd{indx}{:,:};    
    Xarray_centeredd = Xarrayd - repmat(mean(Xarrayd),[size(Xarrayd,1) 1]);
    Xarray_normalizedd = sqrt(glom_movd).*(sqrt(size(Xarrayd,1))*Xarray_centeredd./repmat(sqrt(diag(Xarray_centeredd'*Xarray_centeredd))',[size(Xarray_centeredd,1) 1]));    
    norm_setXd{indx} = setXd{indx};
    norm_setXd{indx}{:,:} = Xarray_normalizedd;
end

merged_normalized = table;
merged_normalized = norm_setXd{1};
for indx = 2:length(setID)
    merged_normalized = [merged_normalized ; norm_setXd{indx}];
end
