%hlid_mds_transform_jackstats_summ: summary plots for hlid_mds_transform_jackstats
%
%runs on results structure from hlid_mds_transform_jackstats
%
%  See also:  HLID_MDS_TRANSFORM_JACKSTATS.
%
rng_state=rng;
if (results.if_frozen~=0) 
    rng('default');
end
colors=rand(results.nstims,3);
rng(rng_state);
%
%magnif factors, as a function of submean, method, dimension and jackknife
%d1: dim-1, d2: first then second magnif then lowest then geomean, d3: full then jackknife, d4: submean d5: method
%
magfacs=zeros(length(results.dimlist)-1,4,1+results.nstims,length(results.submean_use_list),length(results.meth_use_list));
%
for imeth_ptr=1:length(results.meth_use_list)
    imeth=results.meth_use_list(imeth_ptr);
    for isubmean_ptr=1:length(results.submean_use_list)
        isubmean=results.submean_use_list(isubmean_ptr);
        for ijack=0:results.nstims
            if (ijack==0)
                mf_all=results.geo_majaxes{1+isubmean,imeth,1};
            else
                mf_all=results.geo_majaxes_jack_by_stim{1+isubmean,imeth,1,ijack};
            end
            for k=2:results.dimlist(end)
                magfacs(k-1,1:2,ijack+1,1+isubmean,imeth)=mf_all{k,k}.ref.magnifs{1}(1:2)'; %take top two values
                magfacs(k-1,3,ijack+1,1+isubmean,imeth)=mf_all{k,k}.ref.magnifs{1}(end); %lowest value
                magfacs(k-1,4,ijack+1,1+isubmean,imeth)=geomean(mf_all{k,k}.ref.magnifs{1}); %lowest value
            end
        end
    end
end
%
for ifig=1:4
    switch ifig
        case 1
           vplot_all=magfacs(:,1,:,:,:)./magfacs(:,2,:,:,:); %ratio of highest to next-highest 
           vplot_name='highest to next-highest';
           yrange=[1 4];
        case 2
           vplot_all=magfacs(:,1,:,:,:)./magfacs(:,2,:,:,:); %ratio of highest to next-highest 
           vplot_name='highest to next-highest, expanded scale';
           yrange=[1 2];
        case 3
           vplot_all=magfacs(:,1,:,:,:)./magfacs(:,3,:,:,:); %ratio of highest to next-highest 
           vplot_name='highest to lowest';
           yrange=[1 10];
        case 4
            vplot_all=magfacs(:,1,:,:,:)./magfacs(:,4,:,:,:); %ratio of highest to geomean
            vplot_name='highest to geomean';
            yrange=[1 4];
    end
    %
    %each submean and method in a panel
    %
    figure;
    set(gcf,'Position',[50 100 1450 800]);
    set(gcf,'Name',vplot_name);
    set(gcf,'NumberTitle','off');
    for imeth_ptr=1:1+length(results.meth_use_list)
        imeth=results.meth_use_list(1+mod(imeth_ptr-1,length(results.meth_use_list))); %last slot for legend
        for isubmean_ptr=1:length(results.submean_use_list)
            isubmean=results.submean_use_list(isubmean_ptr);
            subplot(length(results.submean_use_list),1+length(results.meth_use_list),imeth_ptr+(isubmean_ptr-1)*(1+length(results.meth_use_list)));
            hold on;
            vplot=reshape(vplot_all(:,1,:,1+isubmean,imeth),length(results.dimlist)-1,1+results.nstims); %d1: dim-1, d2: full then jackknife
            for ijack=1:results.nstims
                hp=plot(results.dimlist(2:end),vplot(:,1+ijack)); %jackknifed ratio
                set(hp,'DisplayName',results.stimulus_names_display{ijack},'LineWidth',2);
                set(hp,'Color',colors(ijack,:));
            end
            hp=plot(results.dimlist(2:end),vplot(:,1),'k');
            set(hp,'DisplayName','full','LineWidth',2);
            hp=plot(results.dimlist(2:end),geomean(vplot(:,2:end),2),'k:');
            set(hp,'DisplayName','geomean(jack)','LineWidth',2);
            %
            set(gca,'XLim',[1.5 results.dimlist(end)+0.5]);
            set(gca,'XTick',results.dimlist(2:end));
            xlabel('dim');
            set(gca,'YLim',yrange);
            set(gca,'YScale','log');
            ylabel('ratio');
            if (imeth_ptr==1+length(results.meth_use_list))
                legend;
                title('legend');
            else
                title(sprintf('sm=%1.0f %s',isubmean,results.meth_names_short{imeth}));
            end
        end %isubmean_ptr
    end %imeth_ptr
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,vplot_name);
    axis off;
    rbase=results.geo{1+results.submean_use_list(1),results.meth_use_list(1)}{1};
    axes('Position',[0.01,0.03,0.01,0.01]); %for text
    text(0,0,sprintf('ref: %s',rbase.ref_file),'Interpreter','none');
    axis off;
    axes('Position',[0.01,0.05,0.01,0.01]); %for text
    text(0,0,sprintf('adj: %s',rbase.adj_file),'Interpreter','none');
    axis off;
    %
    %each dimension and submean in a panel
    %
    figure; %each submean and method in a panel
    set(gcf,'Position',[50 100 1450 800]);
    set(gcf,'Name',vplot_name);
    set(gcf,'NumberTitle','off');
    for idim=results.dimlist(2):results.dimlist(end)
        for isubmean_ptr=1:length(results.submean_use_list)
            isubmean=results.submean_use_list(isubmean_ptr);
            subplot(length(results.dimlist)-1,length(results.submean_use_list),isubmean_ptr+(idim-results.dimlist(2))*length(results.submean_use_list)); ...
            hold on;
            vplot2=reshape(vplot_all(idim-results.dimlist(2)+1,1,:,1+isubmean,:),[1+results.nstims,results.nmeths])';  %d1: method d2: full then jackknife
            for ijack=1:results.nstims
                hp=plot(vplot2(results.meth_use_list,1+ijack)); %jackknifed ratio
                set(hp,'DisplayName',results.stimulus_names_display{ijack},'LineWidth',2);
                set(hp,'Color',colors(ijack,:));
            end
            hp=plot(vplot2(results.meth_use_list,1),'k');
            set(hp,'DisplayName','full','LineWidth',2);
            hp=plot(geomean(vplot2(results.meth_use_list,2:end),2),'k:');
            set(hp,'DisplayName','geomean(jack)','LineWidth',2);
            set(gca,'XLim',[0.5 length(results.meth_use_list)+0.5]);
            set(gca,'XTick',[1:length(results.meth_use_list)]);
            set(gca,'XTickLabel',results.meth_names_short(results.meth_use_list));
            set(gca,'YLim',yrange);
            set(gca,'YScale','log');
            ylabel('ratio');
            title(sprintf('sm=%1.0f dim %1.0f',isubmean,idim));
        end %isubmean_ptr
    end %idim
    axes('Position',[0.01,0.01,0.01,0.01]); %for text
    text(0,0,vplot_name);
    axis off;
    rbase=results.geo{1+results.submean_use_list(1),results.meth_use_list(1)}{1};
    axes('Position',[0.01,0.03,0.01,0.01]); %for text
    text(0,0,sprintf('ref: %s',rbase.ref_file),'Interpreter','none');
    axis off;
    axes('Position',[0.01,0.05,0.01,0.01]); %for text
    text(0,0,sprintf('adj: %s',rbase.adj_file),'Interpreter','none');
    axis off;
end %ifig
