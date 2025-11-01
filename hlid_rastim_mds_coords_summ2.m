%hlid_rastim_mds_coords_summ2: special-purpose script to summarize output
%of several runs of hlid_rastim_mds_coords_demo, looking at 
% fraction of variance explained
%
%   See also:  HLID_RASTIM_MDS_COORDS_DEMO, PSG_KNIT_STATS_PLOT, PSG_ALIGN_STATS_PLOT
%  HLID_RASTIM_MDS_COORDS_SUMM.
%
if ~exist('runs')
    runs=cell(1,4);
    runs{1}.filename='./plots/hlid_rastim_mds_orn-megamat_28Oct25.mat';
    runs{1}.label='orn-megamat';
    runs{1}.color='r';
    runs{2}.filename='./plots/hlid_rastim_mds_kc-megamat_28Oct25.mat';
    runs{2}.label='kc-megamat';
    runs{2}.color='k';
    runs{3}.filename='./plots/hlid_rastim_mds_kc-tnt3c_28Oct25.mat';
    runs{3}.label='kc-tnt3c';
    runs{3}.color='b';
    runs{4}.filename='./plots/hlid_rastim_mds_kc-tntlabel_28Oct25.mat';
    runs{4}.label='kc-tntlabel';
    runs{4}.color='c';    
end
nruns=length(runs);
d=cell(1,nruns);
for irun=1:nruns
    d{irun}=load(runs{irun}.filename);
    disp(sprintf(' data for run %2.0f read from %s',irun,runs{irun}.filename));
end
meth_labels=cell(1,length(d{1}.r.meths));
nmeths=length(meth_labels);
for imeth=1:nmeths
    meth_labels{imeth}=d{1}.r.meths{imeth}.name_short;
end
disp('method labels');
disp(meth_labels);
ra_setup=d{1}.knit_stats_setup;
if ~exist('shuff_quantiles')
    ra_setup.shuff_quantiles=[0.01 0.05 0.5 0.95 0.99];
else
    ra_setup.shuff_quantiles=shuff_quantiles;
end
disp(ra_setup);
%
tstring='rms var unexp/rms var avail, submean: 0 (top) 1 (bottom)';
figure;
set(gcf,'Position',[100 100 1400 800]);
set(gcf,'NumberTitle','off');
set(gcf,'Name',tstring);
if_submean=1;
for submean=0:if_submean
    for imeth=1:nmeths
        hl=cell(0);
        ht=[];
        for irun=1:nruns
            subplot(1+if_submean,nmeths,imeth+submean*nmeths);
            ra=d{irun}.r.knit_stats{imeth,1+submean};
            hp=plot(1:ra_setup.dim_list_in_max,ra.rmsdev_overall(:,1)./ra.rmsavail_overall(:,1),'k');
            set(hp,'Color',runs{irun}.color);
            hl=[hl,hp];
            ht=strvcat(ht,runs{irun}.label);
            hold on;
        end %irun
        set(gca,'XTick',1:ra_setup.dim_list_in_max);
        set(gca,'XLim',[0 ra_setup.dim_list_in_max]);
        xlabel('dim');
        set(gca,'YLim',[0 0.6]);
        ylabel('rms unexp/rms avail');
        title(meth_labels{imeth});
    end %imeth
    legend(hl,ht,'Location','Best','FontSize',7);
end %submean
axes('Position',[0.01,0.04,0.01,0.01]); %for text
text(0,0,tstring,'FontSize',8);
axis off;
