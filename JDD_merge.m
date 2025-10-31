% JDD merge. 
% Getting this working without going through the UI.
% Going to open each of the files in the first two sets.

% dependencies : fileToRaw.m, checkConsist.m

opts = struct;
opts.rep_char = '+';
opts.trial_repeats = 3;
opts.set_names = {'kiTC','meTC','moSK','vaTC'};
opts.kdf = {'max_peak','mean_peak'}; % Known data fields.
opts.suppressoutput = false;
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
[~,Sall] = checkConsist(Sraw,opts);

% First, I want to check that each stimulus has a non-NaN value somewhere
% within a set. If found, that stimulus is removed.
% The presence of the glomeruli across sets is determined. the glomeruli
% used appear in a minimum number of files. 
Strimmed = lookForSetWideHoles(Sall,opts);

% Fill in the remaining holes.
% The calls the afalwt commands
Sfilled = fillInNaNs(Strimmed,opts);

imagesc(Sfilled{1}{1}{:,:})
set(gca,'XTick',[1:length(Sfilled{1}{1}.Properties.VariableNames)]);
set(gca,'XTickLabel',Sfilled{1}{1}.Properties.VariableNames)


% Generate the plots 
%for setindx = 1:numSets
%    numFiles = length(Sfilled{setindx});
%    for ifig=1:3
%        switch ifig
%            case 1
%                figname='all raw data';
%                for fileindx = 1:numFiles
%                    resps_plot(:,:,fileindx) = Sall{setindx}{fileindx}{:,:};
                %glomeruli_plot=[1:nglomeruli];
%            case 2
%                figname='raw data from selected glomeruli';
                %resps_plot=resps_gu;
                %glomeruli_plot=glomeruli_use{iset};
%            case 3
%                figname='raw data from selected glomeruli with missing data filled in';
                %resps_plot=resps_gu_filled;
                %glomeruli_plot=glomeruli_use{iset};
%        end
%        figname=cat(2,sprintf('set %2.0f: ',iset),figname);
%        figure;
%        set(gcf,'Position',[50 100 1800 800]);
%        set(gcf,'NumberTitle','off');
%        set(gcf,'Name',figname);
%        [nr,nc]=nicesubp(nfiles_use(iset));
%        for ifile_ptr=1:nfiles_use(iset)
%            ifile=files_use{iset}(ifile_ptr);
%            subplot(nr,nc,ifile_ptr);
%            imagesc(resps_plot(:,:,ifile_ptr));
%            minmax=[min(min(resps_plot(:,:,ifile_ptr),[],'omitnan')),max(max(resps_plot(:,:,ifile_ptr),[],'omitnan'))];
%            title_string=sprintf('set %1.0f file %2.0f: %s  [%6.3f %6.3f]',iset,ifile,strrep(filenames_short{iset}{ifile},'.mat',''),minmax);
%            title(title_string,'Interpreter','none');
%            set(gca,'FontSize',7);
%            set(gca,'XTick',[1:length(glomeruli_plot)]);
%            set(gca,'XTickLabel',glomeruli(glomeruli_plot));
%            set(gca,'YTick',[1:nstims(iset)]);
%            set(gca,'YTickLabel',stim_labels_set{iset});
%        end
%        axes('Position',[0.01,0.02,0.01,0.01]); %for text
%        text(0,0,figname,'Interpreter','none');
%        axis off;
%    end %ifig








