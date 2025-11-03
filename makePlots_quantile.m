function status = makePlots_quantile(resps_set,resp_range,opts)
% JDD 11-3 
% Make the quantile plots
numSets = length(resps_set);
figure;
set(gcf,'Position',[50 100 800 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name','histogram');


for setindx=1:numSets
    data = resps_set{setindx}{:,:};
    if(~opts.suppressoutput)
        disp(sprintf('quantiles for set %1.0f',setindx))
    end
    subplot(2,numSets,setindx);
    hist(data(:),opts.hist_bins);
    xlabel('response')
    ylabel('counts');
    set(gca,'XLim',resp_range);
    title(sprintf('set %s',opts.set_names{setindx}))
    quantiles=quantile(data(:),opts.hist_quantiles);
    for k=1:length(opts.hist_quantiles)
        axes('Position',[0.01+(setindx-1)/numSets,0.06+0.04*k,0.01,0.01]); %for text
        qt=sprintf(' quantile %5.3f: %7.3f',opts.hist_quantiles(k),quantiles(k));
        if(~opts.suppressoutput)
            disp(qt)
        end
        text(0,0,qt,'Interpreter','none');
        axis off;
    end
end
axes('Position',[0.01,0.02,0.01,0.01]);
text(0,0,sprintf('restore size=%2.0f, subtract mean=%2.0f',opts.restore_size,opts.submean));
axis off;

status = 'Did it plot?';
end

