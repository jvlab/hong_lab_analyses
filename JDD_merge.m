% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

% dependencies : fileToRaw.m, checkConsist.m, lookForSetWideHoles,
% fillInNaNs, makePlots_1, calcResp, makePlots_quantile.

opts = struct;
opts.rep_char = '+';
opts.trial_repeats = 3;
opts.set_names = {'kiTC','meTC','moSK','vaTC'};
opts.kdf = {'max_peak','mean_peak'}; % Known data fields.
opts.suppressoutput = true;
opts.interactive = false;
opts.restore_size = true;
opts.submean = false;
opts.hist_quantiles = [0.05 .25 .5 .75 .96];
opts.hist_bins = 50;

hlid_setup;
% I am assuming that the files that are part of a set are in a separate
% directory. There is not name checking, if a data file is in the folder
% the code will attempt to load it into the set.
%
opts.set_names = {'kiTC','meTC','moSK','vaTC'};
Sraw{1} = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Sraw{2} = fileToRaw('../orn_terminals_Oct25/megamat17');
Sraw{3} = fileToRaw('../orn_terminals_Oct25/monat');
Sraw{4} = fileToRaw('../orn_terminals_Oct25/validation2');







% Switch the if_target arrays to include only the diagnostic odors. 
% I think this is pretty simple to do
%{
for setindx = 1:length(Sraw)
    numFiles = length(Sraw{setindx});
    for fileindx = 1:numFiles
        numTargets = length(Sraw{setindx}{fileindx}.response_amplitude_stim.is_target);
        Sraw{setindx}{fileindx}.response_amplitude_stim.is_target=...
            logical(ones(1,numTargets)-Sraw{setindx}{fileindx}.response_amplitude_stim.is_target);
        numTargets = length(Sraw{setindx}{fileindx}.trial_info.is_target);
        Sraw{setindx}{fileindx}.trial_info.is_target=...
            logical(ones(1,numTargets)-Sraw{setindx}{fileindx}.trial_info.is_target);
    end
end
%}

for setindx = 1:length(Sraw)
    numFiles = length(Sraw{setindx});
    for fileindx = 1:numFiles
        numTargets = length(Sraw{setindx}{fileindx}.response_amplitude_stim.is_target);
        Sraw{setindx}{fileindx}.response_amplitude_stim.is_target=...
            ones(1,numTargets);
        numTargets = length(Sraw{setindx}{fileindx}.trial_info.is_target);
        Sraw{setindx}{fileindx}.trial_info.is_target=...
            ones(1,numTargets);
    end
end



%
%{
opts.set_names = {'kiTC','meTC','vaTC'};
Sraw{1} = fileToRaw('../orn_terminals_Oct25/kiwimix_and_controlmix');
Sraw{2} = fileToRaw('../orn_terminals_Oct25/megamat17');
%Sraw{3} = fileToRaw('../orn_terminals_Oct25/monat');
Sraw{3} = fileToRaw('../orn_terminals_Oct25/validation2');
%}
% Check stimulus consistency across files within a set,
% and glomerus consistency across sets.
[~,Sall] = checkConsist(Sraw,opts);

findRepeatStimuli(Sall);

error('stop')
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
% create the resps_set table. These are the responses, which are taken to
% be the slope of the regression in the afalwt fit.
[resps_set,resp_range] = calcResp(Sfilled,afalwt_fit,opts);

% Generate the quantile plots
makePlots_quantile(resps_set,resp_range,opts);

% Merge the sets into an intersection and union of glomeruli.
merged_data = mergeSets(resps_set);

% All good to here. The merged sets are in the merged data cell array.


for icombm = 1:2
    numSets = length(Sfilled);
    dsid = cell(numSets,1);
    for setindx = 1:numSets
        numFiles = length(Sfilled{setindx});
        for fileindx = 1:numFiles
            dsid_this = Sraw{setindx}{fileindx}.meta.title;
            dsid_this=strrep(dsid_this,'/','-');
            dsid_this=strrep(dsid_this,'\','-');
            dsid_this=strrep(dsid_this,'_','-');
            while contains(dsid_this,'--')
                dsid_this=strrep(dsid_this,'--','-');
            end
        %
            dsid{setindx}=[dsid{setindx};dsid_this];
        end
    end
    f_base=struct;
    for setindx=1:numSets
        f_base.metadata{setindx}=Sraw{setindx}{1}.meta; %original metadata from Hong Lab
        f_base.dsid{setindx}=dsid{setindx}; %data set ID, with special chars turned into -
    end
    f_base.resps=resps_set; %original responses
    if ~opts.restore_size
        f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in'; %original field for responses from Hong Lab
    else
        f_base.coord_opts.resp_type='response_amplitude_stim, after affine filling in and restoring to original size';
    end

    % This is just to get by. I need to figure out why there are four different
    % lists of stimuli.
    stimulusNames = merged_data{icombm}.Properties.RowNames;
    f_base.stimulus_names_orig=char(stimulusNames); %original stimulus names
    f_base.stim_labels_orig=char(stimulusNames); %shortened names for plotting
    f_base.stimulus_names=f_base.stimulus_names_orig;
    f_base.stim_labels=f_base.stim_labels_orig;

    maxdim = min(length(merged_data{icombm}.Properties.RowNames),length(merged_data{icombm}.Properties.VariableNames));
    maxdim_use = maxdim;
    da = merged_data{icombm}{:,:};
    [f,s_diag_all,u_full,v_full,s_full,coords_all]=hlid_coords_svd(f_base,da,maxdim,maxdim_use,opts.submean);

    f_combm{icombm}=f;
    s_diag_all_combm{icombm}=s_diag_all;
    u_full_combm{icombm}=u_full;
    v_full_combm{icombm}=v_full;
    s_full_combm{icombm}=s_full;
    coords_all_combm{icombm}=coords_all;
    %
    f.glomeruli_names=merged_data{icombm}.Properties.VariableNames; %names of glomeruli in each dataset
    f.glomeruli_subset_use=merged_data{icombm}.Properties.VariableNames; %glomeruli used in each subset
    f.glomeruli_combined=merged_data{icombm}.Properties.VariableNames;
    f.roi_names=merged_data{icombm}.Properties.VariableNames;

    for iset = 1:numSets
        resps_combined_set = merged_data{icombm}{resps_set{iset}.Properties.RowNames,:};
        maxdim_set = min(size(resps_combined_set));
        [fset,s_diag_set{iset,icombm},u_full_set{iset,icombm},v_full_set{iset,icombm},s_full_set{iset,icombm}]=...
            hlid_coords_svd(struct(),resps_combined_set,maxdim_set,maxdim_set,opts.submean);    
    end

    resps = merged_data{icombm}{:,:};
    roi_names = merged_data{icombm}.Properties.VariableNames;
    stim_labels = merged_data{icombm}.Properties.RowNames;
    dsid_show='This needs to exist';
    if_submean = opts.submean;
    for setindx = 1:numSets
        nstims(setindx) = length(resps_set{setindx}.Properties.RowNames);
    end
    maxdim_allowed = maxdim;
    if(icombm == 1)
        comb_label = 'intersection';
    else
        comb_label = 'union';
    end

    hlid_coords_plot;
    axes('Position',[0.5,0.02,0.01,0.01]); %for text
    text(0,0,comb_label,'Interpreter','none');
    axis off;

    figname=sprintf('combined sets via %s',comb_label);
    figure;
    set(gcf,'Position',[50 100 1800 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',figname);
    %full array of merged responses
    subplot(1,2,1)
    imagesc(merged_data{icombm}{:,:},resp_range);
    axis equal;
    axis tight;
    hold on;
    for iset=1:numSets-1
        stimct=sum(nstims(1:iset));
        plot([0 length(merged_data{icombm}.Properties.VariableNames)]+[0.5 0.5],stimct+[0.5 0.5],'k','LineWidth',2);
    end
    title_string=figname;
    set(gca,'FontSize',7);
    set(gca,'XTick',[1:length(merged_data{icombm}.Properties.VariableNames)]);
    set(gca,'XTickLabel',merged_data{icombm}.Properties.VariableNames);
    set(gca,'YTick',[1:sum(nstims)]);
    set(gca,'YTickLabel',merged_data{icombm}.Properties.RowNames);
    title(title_string,'Interpreter','none','FontSize',10);
    colorbar;
    %
    %compare PCs in individual sets and full set
    %
    ncols=2*numSets;
    nrows=2;
    for iset=1:numSets
        subplot(nrows,ncols,ncols/2+iset);
        dots=v_full_combm{icombm}'*v_full_set{iset,icombm};
        imagesc(abs(dots),[0 1]);
        xlabel(sprintf('eivec, set %1.0f',iset));
        ylabel('combined eivec');
        axis equal;
        axis tight;
        title('|dots|');
    end
    %
    %compare PCs in individual sets with each other
    %
    nset_pairs=numSets*(numSets-1)/2;
    ncols=2*max(2,nset_pairs);
    ipair=0;
    for jset=2:numSets
        for iset=1:jset-1
            ipair=ipair+1;
            subplot(2,ncols,ncols+ncols/2+ipair);
            dots=v_full_set{jset,icombm}'*v_full_set{iset,icombm};
            imagesc(abs(dots),[0 1]);
            xlabel(sprintf('eivec, set %1.0f',iset));
            ylabel(sprintf('eivec, set %1.0f',jset));
            axis equal;
            axis tight;
            title('|dots|');
        end
    end
    %
    axes('Position',[0.01,0.02,0.01,0.01]); %for text
    text(0,0,figname_raw,'Interpreter','none');
    axis off;
    axes('Position',[0.5,0.02,0.01,0.01]); %for text
    text(0,0,title_string,'Interpreter','none');
    axis off;
end
