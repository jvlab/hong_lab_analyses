function merged_data = mergeSets(S)
% JDD 11-7
% Takes the set response tables and combines them into an intersection and
% a union of glomeruli
numSets = length(S);
%stimulus_list = {};
%for setindx = 1:numSets
%    stimulus_list = [stimulus_list; S{setindx}{1}.Properties.RowNames];
%end
% At this time, there is no need to check for duplicates in this list.

%nglomeruli = length(S{1}{1}.Properties.VariableNames);

glomeruli_combined=cell(2,1);
glomeruli_combined{1}=S{1}.Properties.VariableNames;
glomeruli_combined{2}=[];
for setindx=1:numSets
    glom_use_set = S{setindx}.Properties.VariableNames;
    glomeruli_combined{1}=intersect(glomeruli_combined{1},glom_use_set);
    glomeruli_combined{2}=union(glomeruli_combined{2},glom_use_set);
end

% Intersection is easy
merged_data = cell(2,1);
merged_data{1} = {};
merged_data{2} = {};

for setindx = 1:numSets    
    resps_concat = S{setindx}(:,glomeruli_combined{1});
    merged_data{1} = [merged_data{1}; resps_concat];
end

% Union is harder, but not much.
resps_tmp = S;
for setindx = 1:numSets
    [numStim,~] = size(resps_tmp{setindx});
    gloms_diff = setdiff(glomeruli_combined{2},resps_tmp{setindx}.Properties.VariableNames);
    resps_tmp{setindx}{:,gloms_diff} = zeros(numStim,length(gloms_diff));
    resps_concat = resps_tmp{setindx}(:,glomeruli_combined{2});
    merged_data{2} = [merged_data{2}; resps_concat];    
end

end

