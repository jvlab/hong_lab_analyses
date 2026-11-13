% hlid_vi_coords_knit_heatmaps_rs: create heatmaps from volumetric imaging data and consensus across preps
%
% Based on hlid_vi_coords_knit_auto_rs, but, rather than knitting and saving the consensus files, this creates heatmaps
% of the distances, for the consensus and comparing the individual flies with each other and the consensus
% 
% Retrieves coordinates for individual preps from results_file, created by hlid_vi_pcafilt_coords_auto
% Recreates the consensus calculation, because hlid_vi_write_consensus_auto_rs, which created the consensus from these files, 
% does not save components into coords_knit_auto_rs_avg.mat
%
% Plots heatmaps for:
%   distances in consensus alignment
%   distances between different individual preps, after aligning to consensus (not after aligning to each other)
%   distances between within individual preps, then averaged
%   distances between individual preps and consensus
%
% See hlid_heatmaps_rs for simple heatmap creation
%
%  Several choices for ordering of stimuli (see label_orders)
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_PCAFILT_COORDS_AUTO,
%   HLID_VI_COORDS_KNIT_RS, HLID_METHS_DEFINE,
%   ZHENG_APL_EMBED_PLOT_RS, RS_GET_COORDSETS, RS_KNIT_COORDSETS, RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_XFORM_SPECIFY, RS_XFORM_APPLY
%   HLID_VI_COORDS_KNIT_AUTO_RS, HLID_HEATMAPS_RS.
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
if ~exist('results_file') results_file='hlid_vi_pcafilt_coords_auto_11Sep26.mat'; end
if ~exist('coord_path') coord_path='./data/kc_vi'; end
results_file=getinp('file with results from hlid_vi_pcafilt_coords_auto','s',[],results_file);
coord_path=getinp('path to coordinate files created by hlid_vi_pcafilt_coords_auto','s',[],strrep(coord_path,'\','/'));
coord_path=strrep(strrep(coord_path,'/',filesep),'\',filesep);
%
load(results_file,'results');
results_dims=size(results);
nsets=results_dims(1);
nvariants=prod(results_dims(2:end));
results_res=reshape(results,nsets,nvariants);
disp(sprintf('loaded results from %s, with dimensions',results_file));
%
disp(results_dims)
disp(sprintf('nsets=%3.0f  nvariants available=%5.0f',nsets,nvariants));
disp(sprintf('first coord file of first variant: %s',results_res{1,1}.coord_file));
disp(sprintf(' last coord file of first variant: %s',results_res{end,1}.coord_file));
disp(sprintf('first coord file of  last variant: %s',results_res{1,end}.coord_file));
disp(sprintf(' last coord file of  last variant: %s',results_res{end,end}.coord_file));
%
disp(' ');
for k=1:nsets
    disp(sprintf(' coord file %2.0f of %2.0f: %s',k,nsets,results_res{k,1}.coord_file));
end
%
for k=1:size(results,2)
    sfilts(1,k)=results{1,k,1,1,1,1}.sfilt;
end
disp('spatial filter values');
disp(sfilts);
sfilt_list=getinp('pointer(s) to selection(s)','d',[1 length(sfilts)],1);
%
for k=1:size(results,3)
    resp_measures{1,k}=results{1,1,k,1,1,1}.resp_measure;
end
disp('responses measures');
disp(resp_measures);
resp_measure_list=getinp('pointer(s) to selection(s)','d',[1 length(resp_measures)],1);
%
for k=1:size(results,4)
    pcrits(1,k)=results{1,1,1,k,1,1}.pcrit;
end
disp('pcrit values');
disp(pcrits);
pcrit_list=getinp('pointer(s) to selection(s)','d',[1 length(pcrits)],1);
%
for k=1:size(results,5)
    meth_strings{1,k}=results{1,1,1,1,k,1}.meth_string;
end
disp('methods');
disp(meth_strings);
meth_string_list=getinp('pointer(s) to selection(s)','d',[1 length(meth_strings)],1);
%
for k=1:size(results,6)
    submeans(1,k)=results{1,1,1,1,1,k}.submean;
end
disp('submean');
disp(submeans);
submean_list=getinp('pointer(s) to selection(s)','d',[1 length(submeans)],1);
%
if ~exist('dim_max_in') dim_max_in=5; end
if ~exist('nshuffs') nshuffs=0; end
%
%determine which methods should be knitted with normalized scales
meths=hlid_meths_define;
meth_strings=cell(1,length(meths));
for imeth=1:length(meths)
    meth_strings{imeth}=meths{imeth}.name_file;
end
%assign colors and marker styles
meth_colors=cell(1,length(meths));
meth_colors(contains(meth_strings,'euc_'))=repmat({'k'},1,sum(contains(meth_strings,'euc_')));
meth_colors(contains(meth_strings,'cos_'))=repmat({'c'},1,sum(contains(meth_strings,'cos_')));
meth_colors(contains(meth_strings,'pears_'))=repmat({'m'},1,sum(contains(meth_strings,'pears_')));
%
meth_markers=cell(1,length(meths));
meth_markers(contains(meth_strings,'_svd'))=repmat({'o'},1,sum(contains(meth_strings,'_svd')));
meth_markers(contains(meth_strings,'_mds'))=repmat({'s'},1,sum(contains(meth_strings,'_mds')));
meth_markers(contains(meth_strings,'_sim'))=repmat({'.'},1,sum(contains(meth_strings,'_sim')));
meth_markers(contains(meth_strings,'_chord'))=repmat({'+'},1,sum(contains(meth_strings,'_chord')));
meth_markers(contains(meth_strings,'_ang'))=repmat({'x'},1,sum(contains(meth_strings,'_ang')));
%
if_allow_scales=contains(meth_strings,'euc');
meth_legs=cell(1,length(meths));
if_ok=0;
while (if_ok==0)
    meth_legs=meth_strings;
    for imeth=1:length(meths)
        if if_allow_scales(imeth)
            meth_legs{imeth}=cat(2,meth_legs{imeth},' scaling');
        end
        disp(sprintf(' method %2.0f is %20s knitted with allow_scale set to %1.0f;  legend is  %s',imeth,meth_strings{imeth},if_allow_scales(imeth),meth_legs{imeth}));
    end
    if_ok=getinp('1 if ok','d',[0 1],1);
    if ~if_ok
        if_allow_scales=getinp('new values','d',[0 1],if_allow_scales);
    end
end
%
dim_max_in=getinp('maximum dimension to consider for consensus','d',[1 24],dim_max_in);
dim_list=getinp('dimension(s) to plot in heatmap','d',[1 dim_max_in],[2 3]);
nshuffs=0; %no statistics
if_remove_pipeline=getinp('1 to remove pipeline field from results_knit consensus data','d',[0 1],1);
%
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.if_pca=0; %if_c2p handled later
opts_knit.max_niters=1000;
opts_knit.pcon_init_method=0;
opts_knit.if_stats=0;
opts_knit.if_plot=0;
opts_knit.dim_max_in=dim_max_in;
opts_knit.nshuffs=nshuffs;
opts_knit.if_frozen=1;
opts_knit.if_log=0;
%
results_knit=cell(1,nvariants);
warnings=[];
%
%stimulus names will be reordered after alignment
raw_stim_labels=results{1,1,1,1,1,1}.stim_labels; %raw file stim labels, in order of coordinates
raw_stim_names=results{1,1,1,1,1,1}.stims.names;
raw_stim_names_short=results{1,1,1,1,1,1}.stims.names_short;
display_order=results{1,1,1,1,1,1}.stims.display_order; %pointers into raw_stim_labels in order of pc and timecourse display of hlid_vi_viewpcs
%
label_orders=struct;
label_orders.vi_behavioral=raw_stim_names_short(display_order); %stimulus labels in the order to be plotted in the heatmap
label_orders.vi_tca={'CO2','p-cre','menth','2h','geos','1,4d','3-8ol','IaA','ma','pfo','mint','ms','H2O','3mt1p','t2h','2,5d','1-6ol','EtOH','4-mch','ban','ga','aa','wine','ACV'};
label_orders.alphabetical=sort(raw_stim_names_short);
label_order_list=fieldnames(label_orders);
for k=1:length(label_order_list)
    disp(sprintf('%1.0f->%s',k,label_order_list{k}));
end
label_order_choice=getinp('choice','d',[1 length(label_order_list)],1);
label_order=label_orders.(label_order_list{label_order_choice});
%
nstims=length(raw_stim_labels);
consensus_ticks=[1:nstims];
indiv_ticks=(1+nsets)/2+[0:nsets:(nstims-1)*nsets];
%
nplots=4;
[nr_plot,nc_plot]=nicesubp(nplots,0.7);
%
for sfilt_ptr=1:length(sfilt_list)
    sfilt=sfilt_list(sfilt_ptr);
    for resp_measure_ptr=1:length(resp_measure_list)
        resp_measure=resp_measure_list(resp_measure_ptr);
        for pcrit_ptr=1:length(pcrit_list)
            pcrit=pcrit_list(pcrit_ptr);
            for method_ptr=1:length(meth_string_list)
                method=meth_string_list(method_ptr);
                for submean_ptr=1:length(submean_list)
                    submean=submean_list(submean_ptr);
                    var_string=sprintf(' spat filt: %1.0f; pcrit: %5.4f; resp measure %10s; dimred method %12s; submean: %1.0f',...
                        sfilt-1,pcrits(pcrit),resp_measures{resp_measure},meth_strings{method},submeans(submean));
                    disp(var_string);
                    fullnames=cell(1,nsets);
                    for iset=1:nsets
                        r=results{iset,sfilt,resp_measure,pcrit,method,submean};
                        fullnames{iset}=cat(2,coord_path,filesep,r.coord_file);
                    end
                    aux=struct;
                    aux.nsets=nsets;
                    [data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
                    disp(sprintf('first file %s',fullnames{1}));
                    disp(sprintf(' last file %s',fullnames{nsets}));
                    %align and knit
                    [data_aligned,aux_align_out]=rs_align_coordsets(data_read,setfields(struct(),{'opts_align','opts_check'},{opts_align,opts_check}));
                    disp('aligned');
                    %
                    meth_string=r.meth_string;
                    imeth=find(contains(meth_strings,meth_string)>0);
                    if length(imeth)==1
                        opts_knit.allow_scale=if_allow_scales(imeth);
                        aux.opts_knit=opts_knit;
                        [data_consensus,aux_knit_out]=rs_knit_coordsets(data_aligned,setfields(struct(),{'opts_knit','opts_check'},{opts_knit,opts_check}));
                        disp(sprintf('knit with dim_max_in=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',opts_knit.dim_max_in,opts_knit.pcon_init_method,opts_knit.allow_scale));
                        if if_remove_pipeline
                            data_consensus.sets{1}=struct();
                        end
                    else
                        wmsg=sprintf('method %s not recognized',meth_string);
                        warning(wmsg);
                        warnings=strvcat(warnings,wmsg);
                    end
                    dptrs=zeros(1,nstims);
                    %determine display order (should not depend on data, but needs to be done after alignment and knitting)
                    for istim=1:nstims
                        stim_match=strmatch(label_order{istim},data_consensus.sas{1}.typenames,'exact');
                        if length(stim_match)==1
                            dptrs(istim)=stim_match;
                        else
                            error(sprintf('stimulus %s not found',label_order{istim}));
                        end
                    end
                    for idim=dim_list
                        tstring=cat(2,sprintf('dim %1.0f: ',idim),var_string);
                        figure;
                        set(gcf,'Position',[50 50 1200 800]);
                        set(gcf,'NumberTitle','off');
                        set(gcf,'Name',tstring);
                        for plot_type=1:nplots
                            switch plot_type
                                case 1 %distances in the consensus
                                    dists=sqrt(cootodsq(data_consensus.ds{1}{idim}));
                                    dists_plot=dists(dptrs,dptrs);
                                    dist_title='consensus vs consensus';
                                    xticks=consensus_ticks;
                                    yticks=consensus_ticks;
                                case 2 %individual vs individual
                                    dists_ivsi=zeros(nsets,nstims,nsets,nstims);
                                    for iset1=1:nsets
                                        for iset2=1:iset1
                                            dists_ivsi(iset1,:,iset2,:)=reshape(sqrt(cootodsq(aux_knit_out.components.ds{iset1}{idim},aux_knit_out.components.ds{iset2}{idim})),[1 nstims 1 nstims]);
                                            if iset2~=iset1
                                                dists_ivsi(iset2,:,iset1,:)=dists_ivsi(iset1,:,iset2,:);
                                            end
                                        end
                                    end
                                    dists_plot=reshape(dists_ivsi(:,dptrs,:,dptrs),nstims*nsets,nstims*nsets);
                                    dist_title='indiv vs indiv';
                                    xticks=indiv_ticks;
                                    yticks=indiv_ticks;
                                case 3 % indiv vs indiv avg
                                    dists_ivsi_avg=zeros(nstims,nstims,nsets);
                                    for iset=1:nsets
                                        dists_ivsi_avg(:,:,iset)=reshape(dists_ivsi(iset,:,iset,:),[nstims,nstims]);                                      
                                    end
                                    dists_plot=mean(dists_ivsi_avg(dptrs,dptrs,:),3);
                                    dist_title='avg(indiv vs indiv)';
                                    xticks=consensus_ticks;
                                    yticks=consensus_ticks;
                                case 4 % indiv vs consensus
                                    dists_ivsc=zeros(nsets,nstims,nstims);
                                    for iset=1:nsets
                                        dists_ivsc(iset,:,:)=reshape(sqrt(cootodsq(data_consensus.ds{1}{idim},aux_knit_out.components.ds{iset}{idim})),[1 nstims nstims]);
                                    end
                                    dists_plot=reshape(dists_ivsc(:,dptrs,dptrs),nstims*nsets,nstims);
                                    dist_title='indiv vs consensus';
                                    xticks=indiv_ticks;
                                    yticks=consensus_ticks;
                            end
                            subplot(nr_plot,nc_plot,plot_type);
                            imagesc(dists_plot');
                            title(dist_title);
                            set(gca,'XTick',xticks);
                            set(gca,'XTickLabel',label_order); %properly ordered
                            set(gca,'YTick',yticks);
                            set(gca,'YTickLabel',get(gca,'XTickLabel'));
                            axis square;
                            axis tight;
                            colorbar;
                        end %plot_type
                        %
                        axes('Position',[0.01,0.01,0.01,0.01]);
                        text(0,0,tstring,'Interpreter','none');
                        axis off
                        axes('Position',[0.01,0.03,0.01,0.01]);
                        text(0,0,sprintf('label order: %s',label_order_list{label_order_choice}),'Interpreter','none');
                        axis off
                        axes('Position',[0.01,0.05,0.01,0.01]);
                        text(0,0,sprintf('results file: %s',results_file),'Interpreter','none');
                        axis off
                        disp(sprintf('plotted %s',tstring));
                    end %dimension
                end %submean
            end %method
        end %pcrit
    end %resp_measeure
end %sfilt
