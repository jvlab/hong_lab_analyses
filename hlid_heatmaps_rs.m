% hlid_heatmaps_rs: create heatmaps from coordinate files
%
% A streamlined version of hlid_vi_coords_knit_heatmaps_rs.
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_COORDS_KNIT_HEATMAPS_RS.
%
hlid_setup;
%
if ~exist('prefix_remove') prefix_remove='hlid_vi_'; end
%
opts_read=struct();
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.if_log=0;
%
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
opts_align.if_log=0;
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
%
aux=struct;
aux.nsets=0;
fullnames=[];
[data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
nsets=length(data_read.ds);
fullnames=cell(1,nsets);
filenames=cell(1,nsets);
for iset=1:nsets
    filename=data_read.sets{1}.label;
    fullnames{iset}=filename;
    seps=union(find(filename=='/'),find(filename=='\'));
    if ~isempty(seps)
        filename=filename(max(seps)+1:end);
    end
    filenames{iset}=filename;
    disp(sprintf('file %1.0f: %s',iset,filenames{iset}));
end
%
if ~exist('label_orders')
    label_orders=struct;
    label_orders=display_orders; %standard orders for kc clust, kctnt, and kc merge
        %kcclust: {'2h' 'IaA' 'pa' '2-but' 'eb' 'ep' 'aa' 'va' 'B-cit' 'Lin' '6al' 't2h' '1-8ol' '1-5ol' '1-6ol' 'benz' 'ms'}
        %kctnt: {'1-5ol' '1-6ol' '1-8ol' '2-but' '2h' '6al' 'IaA' 'PAA' 'eb' 'ep' 'ms' 'pa' 't2h' 'B-myr' 'euc' '-aPine' 'pfo'}
        %kcmerge: {'2h' 'IaA' 'pa' '2-but' 'eb' 'ep' '6al' 't2h' '1-8ol' '1-5ol' '1-6ol' 'PAA' 'ms' 'B-myr' 'euc' '-aPine' 'pfo'}
        %orders for volumetric imaging
    label_orders.vi_behavioral={'CO2' 'p-cre' 'menth' 'IaA' '1,4d' 'mint' '3-8ol' '1-6ol' 'pfo' 'ms' '2h' '2,5d' 'H2O' 'geos' 'ma' 't2h' 'ga' '3mt1p' 'aa' '4-mch' 'EtOH' 'ban' 'wine' 'ACV'};
    label_orders.vi_tca={'CO2' 'p-cre' 'menth' '2h' 'geos' '1,4d' '3-8ol' 'IaA' 'ma' 'pfo' 'mint' 'ms' 'H2O' '3mt1p' 't2h' '2,5d' '1-6ol' 'EtOH' '4-mch' 'ban' 'ga' 'aa' 'wine' 'ACV'};
    label_orders.vi_alphabetical={'1,4d' '1-6ol' '2,5d' '2h' '3-8ol' '3mt1p' '4-mch' 'ACV' 'CO2' 'EtOH' 'H2O' 'IaA' 'aa' 'ban' 'ga' 'geos' 'ma' 'menth' 'mint' 'ms' 'p-cre' 'pfo' 't2h' 'wine'};
end
label_order_list=fieldnames(label_orders);
for k=1:length(label_order_list)
    labels=label_orders.(label_order_list{k});
    disp(sprintf('%2.0f->%20s (%4.0f labels) : %6s %6s %6s %6s ... %6s %6s %6s %6s',k,label_order_list{k},length(labels),labels{[1 2 3 4 end-3 end-2 end-1 end]}));
end
label_order_choice=getinp('choice (0: native order in file','d',[0 length(label_order_list)],1);
%
dim_list=getinp('dimension(s) to plot in heatmap','d',[1 length(data_read.ds{1})],[2 3]);
%
for iset=1:nsets
    disp(sprintf('processing %s',filenames{iset}));
    typenames=data_read.sas{iset}.typenames;
    if label_order_choice==0
        nstims=length(typenames);
        dptrs=[1:nstims];
        label_order=typenames;
    else
        label_order=label_orders.(label_order_list{label_order_choice});
        nstims=length(label_order);
        dptrs=zeros(1,nstims);
        %determine display order
        for istim=1:nstims
            stim_match=strmatch(label_order{istim},data_read.sas{1}.typenames,'exact');
            if length(stim_match)==1
                dptrs(istim)=stim_match;
            else
                disp(sprintf('no stimulus with label %s is found',labels_order{istim}));
            end
        end %istim
    end
    disp(sprintf(' stimuli in file: %3.0f, to be plotted: %3.0f, in label list but not in file: %3.0f, in file but not in label list: %3.0f',...
        length(typenames),sum(dptrs>0),sum(dptrs==0),length(setdiff([1:length(typenames)],dptrs))));
    for idim=dim_list
        dists=sqrt(cootodsq(data_read.ds{iset}{idim}));
        dists_plot=zeros(nstims,nstims);
        dptrs_found=dptrs(dptrs>0);
        dists_plot(dptrs>0,dptrs>0)=dists(dptrs(dptrs>0),dptrs(dptrs>0));
        %
        tstring=cat(2,sprintf('dim %1.0f: ',idim),filenames{iset});
        figure;
        set(gcf,'Position',[50 50 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        imagesc(dists_plot);
        set(gca,'XTick',[1:nstims]);
        set(gca,'XTickLabel',label_order); %properly ordered
        set(gca,'YTick',[1:nstims]);
        set(gca,'YTickLabel',get(gca,'XTickLabel'));
        axis square;
        axis tight;
        colorbar;
        %
        axes('Position',[0.01,0.01,0.01,0.01]);
        text(0,0,tstring,'Interpreter','none');
        axis off
    end %idim
end
