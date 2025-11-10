%hlid_rastim_mds_coords_summ2s: special-purpose script to summarize output
%of several runs of hlid_rastim_mds_coords_demo, looking at 
% fraction of variance explained
%
% in contrast to hlid_rstim_mds_coords_summ2, this looks at subsample analysis,
%   from hlid_rastim_mds_coords_subsamp.
%
%   See also:  HLID_RASTIM_MDS_COORDS_DEMO, HLID_RASTIM_MDS_COORDS_SUBSAMP,
%  HLID_RASTIM_MDS_COORDS_SUMM2.
%
if ~exist('runs') %files from hlid_rastim_mds_coords_subsamp
    runs=cell(1,4);
    runs{1}.filename='./plots/hlid_rastim_mds_orn-megamat_07Nov25.mat';
    runs{1}.filename_subsamp='./plots/hlid_rastim_mds_subsamp_orn-megamat_09Nov25.mat';
    runs{1}.label='orn-megamat';
    runs{1}.color='r';
    runs{2}.filename='./plots/hlid_rastim_mds_kc-megamat_07Nov25.mat';
    runs{2}.filename_subsamp='./plots/hlid_rastim_mds_subsamp_kc-megamat_09Nov25.mat';
    runs{2}.label='kc-megamat';
    runs{2}.color='k';
    runs{3}.filename='./plots/hlid_rastim_mds_kc-tnt3c_07Nov25.mat';
    runs{3}.filename_subsamp='./plots/hlid_rastim_mds_subsamp_kc-tnt3c_09Nov25.mat';
    runs{3}.label='kc-tnt3c';
    runs{3}.color='b';
    runs{4}.filename='./plots/hlid_rastim_mds_kc-tntlabel_07Nov25.mat';
    runs{4}.filename_subsamp='./plots/hlid_rastim_mds_subsamp_kc-tntlabel_09Nov25.mat';
    runs{4}.label='kc-tntlabel';
    runs{4}.color='c';    
end
nruns=length(runs);
d=cell(1,nruns);
for irun=1:nruns
    d{irun}=load(runs{irun}.filename); %load basic info
    d{irun}.r_subsamp=getfield(load(runs{irun}.filename_subsamp),'r_subsamp'); %load subsamp info
    disp(sprintf(' data for run %2.0f read from %s',irun,runs{irun}.filename_subsamp));
end
meth_labels=cell(1,length(d{1}.r.meths));
nmeths=length(meth_labels);
for imeth=1:nmeths
    meth_labels{imeth}=d{1}.r.meths{imeth}.name_short;
end
disp('method labels');
disp(meth_labels);
ra_setup=d{1}.knit_stats_setup;
%disp(ra_setup);
disp(sprintf('nstims from knit_stats_setup: %2.0f',ra_setup.nstims));
disp(sprintf('nsets: %2.0f',ra_setup.nsets));
%
disp('from r_subsamp:');
dim_max_in=d{1}.r_subsamp.dim_max_in;
subsamp_sizes=d{1}.r_subsamp.subsamp_sizes;
nsizes=d{1}.r_subsamp.nsizes;
for isize=1:nsizes
    disp(sprintf(' subsample into %2.0f files',subsamp_sizes(isize)));
    for irun=1:nruns
        disp(sprintf(' %15s (%2.0f files available): number of subsamples: %5.0f',...
            runs{irun}.label,d{irun}.r_subsamp.nfiles,size(d{irun}.r_subsamp.subsamp_lists{isize},1)));
    end
end
%
for isize=1:nsizes
    subsamp_size=subsamp_sizes(isize);
    tstring=sprintf('subsamp: %2.0f files: rms var unexp/rms var avail',subsamp_size);
    tstring2=cat(2,tstring,' submean: 0 (top) 1 (bottom)');
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
                %plot values without subsamples
                ra=d{irun}.r.knit_stats{imeth,1+submean}; %no subsample
                hp=plot(1:ra_setup.dim_list_in_max,ra.rmsdev_overall(:,1)./ra.rmsavail_overall(:,1),'k:','LineWidth',2);
                set(hp,'Color',runs{irun}.color);
                hold on;
                hl=[hl,hp];
                ht=strvcat(ht,cat(2,runs{irun}.label,', all'));
            end
            for irun=1:nruns
                %plot mean of values with subsamples
                rs=d{irun}.r_subsamp.stats{imeth,1+submean,isize};
                if ~isempty(rs)
                    hp=plot(1:dim_max_in,mean(rs.rmsdev_overall(:,1,:)./rs.rmsavail_overall(:,1,:),3,'omitnan'),'k','LineWidth',1);
                    set(hp,'Color',runs{irun}.color);
                    hl=[hl,hp];
                    ht=strvcat(ht,sprintf('subsampled (%1.0f)',subsamp_size));
                end
            end %irun
            set(gca,'XTick',1:dim_max_in);
            set(gca,'XLim',[0 dim_max_in]);
            xlabel('dim');
            set(gca,'YLim',[0 0.6]);
            ylabel('rms unexp/rms avail');
            title(meth_labels{imeth});
        end %imeth
        legend(hl,ht,'Location','Best','FontSize',7);
    end %submean
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,tstring2,'FontSize',8);
    axis off;
end %isize
