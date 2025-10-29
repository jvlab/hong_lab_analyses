%hlid_rastim_mds_coords_summ: special-purpose script to summarize output
%of several runs of hlid_rastim_mds_coords_demo, looking at consistency across shuffles
%
%   See also:  HLID_RASTIM_MDS_COORDS_DEMO, PSG_KNIT_STATS_PLOT, PSG_ALIGN_STATS_PLOT.
%
if ~exist('runs')
    runs=cell(1,4);
    runs{1}.filename='./plots/hlid_rastim_mds_orn-megamat_28Oct25.mat';
    runs{1}.label='orn-megamat';
    runs{2}.filename='./plots/hlid_rastim_mds_kc-megamat_28Oct25.mat';
    runs{2}.label='kc-megamat';
    runs{3}.filename='./plots/hlid_rastim_mds_kc-tnt3c_28Oct25.mat';
    runs{3}.label='kc-tnt3c';
    runs{4}.filename='./plots/hlid_rastim_mds_kc-tntlabel_28Oct25.mat';
    runs{4}.label='kc-tntlabel';
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
nquantiles=length(ra_setup.shuff_quantiles);
%
if_submean=1;
for iue=1:2 %explained or unexplained rms
    if (iue==1)
        var_string='unexplained'; 
        iue_sign=1; %logic to add unexplained variance
        iue_mult=0; %and not include available variance
    else
        var_string='explained';
        iue_sign=-1; %logic to subtract unexplained variance
        iue_mult=1; %and add explained variance
    end
    for submean=0:if_submean
        tstring=sprintf('%s var, preproc submean=%1.0f',var_string,submean);
        figure;
        set(gcf,'Position',[100 100 1400 900]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        ymax=-Inf;
        for irun=1:nruns
            for imeth=1:nmeths
                subplot(nruns,nmeths,imeth+(irun-1)*nmeths);
                ra=d{irun}.r.knit_stats{imeth,1+submean};
                hl=cell(0);
                hp=plot(1:ra_setup.dim_list_in_max,iue_mult*ra.rmsavail_overall(:,1)+iue_sign*ra.rmsdev_overall(:,1),'k');
                hl=[hl,hp];
                ht='consensus, data';
                hold on;
                if (iue==2)
                    hp=plot(1:ra_setup.dim_list_in_max,ra.rmsavail_overall(:,1),'b');
                    hl=[hl,hp];
                    ht=strvcat(ht,'avail');
                end
                for iq=1:nquantiles
                    switch sign(ra_setup.shuff_quantiles(iq)-0.5)
                        case -1
                            linetype=':';
                        case 0
                            linetype='';
                        case 1
                            linetype='--';
                    end
                    hp_last=plot(1:ra_setup.dim_list_in_max,iue_mult*ra.rmsavail_overall(:,1)+...
                        iue_sign*quantile(ra.rmsdev_overall_shuff(:,1,1,:,1),ra_setup.shuff_quantiles(iq),4),cat(2,'r',linetype));
                    hp_all=plot(1:ra_setup.dim_list_in_max,iue_mult*ra.rmsavail_overall(:,1)+...
                        iue_sign*quantile(ra.rmsdev_overall_shuff(:,1,1,:,2),ra_setup.shuff_quantiles(iq),4),cat(2,'m',linetype));
                    if iq==round((1+nquantiles)/2)
                        hl=[hl,hp_last,hp_all];
                        ht=strvcat(ht,'cons, last','cons, all');
                     end
                end
                set(gca,'XTick',1:ra_setup.dim_list_in_max);
                set(gca,'XLim',[0 ra_setup.dim_list_in_max]);
                xlabel('dim');
                title(meth_labels{imeth});
                ylabel(runs{irun}.label);
                if (irun==nruns) & (imeth==nmeths)
                    legend(hl,ht,'Location','Best','FontSize',7);
                end
                ymax=max(ymax,max(get(gca,'YLim')));
            end %imeth
        end %irun
        %equalize ordinate scale
        for irun=1:nruns
            for imeth=1:nmeths
                subplot(nruns,nmeths,imeth+(irun-1)*nmeths);
                set(gca,'YLim',[0 ymax]);
            end
        end
        axes('Position',[0.01,0.04,0.01,0.01]); %for text
        text(0,0,cat(2,sprintf('%5s, quantiles from %5.0f shuffles: ',tstring,ra_setup.nshuffs),sprintf('%6.4f ',ra_setup.shuff_quantiles)),...
            'FontSize',8);
    axis off;
    end%submean
end %iue
