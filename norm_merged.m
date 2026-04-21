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

norm_setX = cell(length(setID),1);
for indx = 1:length(setX)
    Xarray = setX{indx}{:,:};
    diag(Xarray'*Xarray)'
    Xarray_centered = Xarray - repmat(mean(Xarray),[size(Xarray,1) 1]);
    Xarray_normalized = sqrt(size(Xarray,1))*Xarray_centered./repmat(sqrt(diag(Xarray_centered'*Xarray_centered))',[size(Xarray_centered,1) 1]);
    norm_setX{indx} = setX{indx};
    norm_setX{indx}{:,:} = Xarray_normalized;
end


merged_normalize = table;
merged_normalize = norm_setX{1};
for indx = 2:length(setID)
    merged_normalize = [merged_normalize ; norm_setX{indx}];
end