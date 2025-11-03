function status = makePlots_1(S1,S2,S3)
% JDD 11-3
% Generates first set of plots in hlid_orn_merge2

numSets = length(S3);

for setindx = 1:numSets
    numFiles = length(S3{setindx});
    for ifig=1:3
        switch ifig
            case 1
                figname = 'all raw data';
                S_plot = S1;
            case 2
                figname = 'raw data from selected glomeruli';
                S_plot = S2;
            case 3
                figname = 'raw data from selected glomeruli with missing dat filled in';
                S_plot = S3;
        end
        [numStim,numGlom] = size(S_plot{setindx}{1});
        
        figname=cat(2,sprintf('set %2.0f: ',setindx),figname);
        
        figure;
        set(gcf,'Position',[50 100 1800 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figname);
        [nr,nc]=nicesubp(numFiles);
        
        for ifile_ptr=1:numFiles
            filename = S_plot{setindx}{ifile_ptr}.Properties.Description;
            filename = split(filename,'/');
            filename = filename{end};
            subplot(nr,nc,ifile_ptr);
            imagesc(S_plot{setindx}{ifile_ptr}{:,:});
            minmax=[min(min(S_plot{setindx}{ifile_ptr}{:,:},[],'omitnan')),max(max(S_plot{setindx}{ifile_ptr}{:,:},[],'omitnan'))];
            title_string=sprintf('set %1.0f file %2.0f: %s  [%6.3f %6.3f]',setindx,ifile_ptr,filename,minmax);
            title(title_string,'Interpreter','none');
            set(gca,'FontSize',7);
            set(gca,'XTick',[1:numGlom]);
            set(gca,'XTickLabel',S_plot{setindx}{ifile_ptr}.Properties.VariableNames);
            set(gca,'YTick',[1:numStim]);
            set(gca,'YTickLabel',S_plot{setindx}{ifile_ptr}.Properties.RowNames);
        end
    end
end

status = string('I had to put something there');

end