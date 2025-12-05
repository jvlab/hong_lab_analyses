opts.stimnametype='long';

stimulus_all = sa.typenames;
% Should sort these.

groupsize = [];
groupList = cell(100,1); % Bigger than we need on purpose.
stim_short = stimulus_all;

for stimindx = 1:length(stimulus_all)
    stim_short{stimindx} = stimulus_all{stimindx}(1:end-4);% This removes the 4 character set code
    if stim_short{stimindx}(1:5)=='diag_'
        stim_short{stimindx} = stim_short{stimindx}(6:end);% This removes the diagnostic designation
    end
end

short_stim_names = stim_short;
while length(stim_short)>1    
    stim1 = stim_short{1};
    findx = strcmp(stim_short,stim1);
    if(sum(findx)>1)
        groupsize = [groupsize sum(findx)];
        tmpList = stimulus_all(findx);
        groupList{length(groupsize)} = tmpList;
    end
    stimulus_all(findx) = [];
    stim_short(findx) = [];    
end

groupList(length(groupsize)+1:end)=[];

% Find the first square that is greater than group num
numGroups = length(groupsize);
it = 1;
while 1
    if(it^3<numGroups)
        it = it + 1;
    else 
        break;
    end
end

step_size = it^3/numGroups;

numColors = numGroups;

bracket = cumsum(groupsize);

sasub = cell(numGroups,1);
%{
for group=1:numGroups
    if(group == 1)
        begin_it = 1;
    else
        begin_it = bracket(group-1)+1;
    end
    end_it = bracket(group);
    sasub{group} = sa;
    sasub{group}.typenames = sa.typenames(begin_it:end_it);
    if(group == 1)
        ou=psg_plotcoords(d{3}(begin_it:end_it,:),[1 2 3],sasub{group},[],setfields(opts_plot,{'color_norays','axis_handle'},{myColormap(group,:),[]}));
    else
        ou=psg_plotcoords(d{3}(begin_it:end_it,:),[1 2 3],sasub{group},[],setfields(opts_plot,{'color_norays','axis_handle'},{myColormap(group,:),ou.axis_handle}));
    end
end
%}
groupsToPlot = 1:numGroups;
myColormap = zeros(length(groupsToPlot),3);

myColormap(:,1) = 0.5*cos(pi*[0:length(groupsToPlot)-1]/length(groupsToPlot))+0.5;
myColormap(:,2) = 0.5*sin(2*pi*[0:length(groupsToPlot)-1]/length(groupsToPlot))+0.5;
myColormap(:,3) = 0.5*cos(4*pi*[0:length(groupsToPlot)-1]/length(groupsToPlot))+0.5;

%numClusters = max(idx);
%myColormap = zeros(numClusters,3);
%myColormap(:,1) = linspace(0,1,numClusters);
%myColormap(:,2) = linspace(1,0,numClusters);
%myColormap(:,1) = 0.5*cos(pi*[0:numClusters-1]/numClusters)+0.5;
%myColormap(:,2) = 0.5*sin(2*pi*[0:numClusters-1]/numClusters)+0.5;
%myColormap(:,3) = 0.5*cos(4*pi*[0:numClusters-1]/numClusters)+0.5;


for group=1:numGroups
    sasub{group} = sa;
    if(group == 1)
        begin_it = 1;
    else
        begin_it = bracket(group-1)+1;
    end
    end_it = bracket(group);
    if(opts.stimnametype == 'long')
        sasub{group}.typenames = sa.typenames(begin_it:end_it);
    elseif(opts.stimnametype == 'short')
        sasub{group}.typenames = short_stim_names(begin_it:end_it);
    else
        sasub{group}.typenames = sa.typenames(begin_it:end_it);
    end
end

%color the plot using kmeans clustering
% need a stand alone data set



for group = 1:length(groupsToPlot)
    if(groupsToPlot(group) == 1)
        begin_it = 1;
    else
        begin_it = bracket(groupsToPlot(group)-1)+1;
    end
    end_it = bracket(groupsToPlot(group));
    
    if(group == 1)
        ou=psg_plotcoords(d{3}(begin_it:end_it,:),[1 2 3],sasub{groupsToPlot(group)},[],setfields(opts_plot,{'color_norays','axis_handle'},{myColormap(group,:),[]}));
        %ou=psg_plotcoords(d{3}(begin_it:end_it,:),[1 2 3],sasub{groupsToPlot(group)},[],setfields(opts_plot,{'color_norays','axis_handle'},{myColormap(idx(groupsToPlot(group)),:),[]}));

    else
        ou=psg_plotcoords(d{3}(begin_it:end_it,:),[1 2 3],sasub{groupsToPlot(group)},[],setfields(opts_plot,{'color_norays','axis_handle'},{myColormap(group,:),ou.axis_handle}));
        
        %ou=psg_plotcoords(d{3}(begin_it:end_it,:),[1 2 3],sasub{groupsToPlot(group)},[],setfields(opts_plot,{'color_norays','axis_handle'},{myColormap(idx(groupsToPlot(group)),:),ou.axis_handle}));
        
    end
end




    