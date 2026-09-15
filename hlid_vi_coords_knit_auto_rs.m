% hlid_vi_coords_knit_auto_rs: read volumetric imaging coordinates set files,
% automated knit and compare, based on outputs of hlid_vi_pcafilt_coords_auto.
%
% Uses results from hlid_vi_pcafilt_coords_auto to determine names of coordinate files and the options usef for conversion to coordinates
% Reads these coordinate files
% Consensus information kept in results_knit
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_PCAFILT_COORDS_AUTO,
%   HLID_VI_COORDS_KNIT_RS, HLID_METHS_DEFINE,
%   ZHENG_APL_EMBED_PLOT_RS, RS_GET_COORDSETS, RS_KNIT_COORDSETS, RS_DISP_COORDSETS, RS_CONCAT_COORDSETS, RS_XFORM_SPECIFY, RS_XFORM_APPLY.
%
hlid_setup;
%
if ~exist('axis_view') axis_view=[-37.5000   30.0000]; end
if ~exist('markersize_consensus') markersize_consensus=24; end
if ~exist('markersize_component') markersize_component=16; end
if ~exist('linewidth') linewidth=2; end
if ~exist('callout_amount') callout_amount=0.5; end
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
coord_path=getinp('path to coordinate files created by hlid_vi_pcafilt_coords_auto','s',[],coord_path);
coord_path=strrep(strrep(coord_path,'/',filesep),'\',filesep);
%
load(results_file,'results');
results_dims=size(results);
nsets=results_dims(1);
nvariants=prod(results_dims(2:end));
results_res=reshape(results,nsets,nvariants);
disp('loaded results, with dimensions')
disp(results_dims)
disp(sprintf('nsets=%3.0f  nvariants=%5.0f',nsets,nvariants));
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
for k=1:size(results,4)
    pcrits(1,k)=results{1,1,1,k,1,1}.pcrit;
end
disp('spatial filter values');
disp(sfilts);
disp('pcrit values');
disp(pcrits);
%
if ~exist('dim_max_in') dim_max_in=10; end
if ~exist('nshuffs') nshuffs=100; end
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
    if_ok=getinp('1 if ok','d',[0 1]);
    if ~if_ok
        if_allow_scales=getinp('new values','d',[0 1],if_allow_scales);
    end
end
%
dim_max_in=getinp('maximum dimension to consider','d',[1 24],dim_max_in);
nshuffs=getinp('number of shuffles','d',[0 1000],nshuffs);
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
opts_knit.if_stats=1;
opts_knit.if_plot=0;
opts_knit.dim_max_in=dim_max_in;
opts_knit.nshuffs=nshuffs;
opts_knit.if_frozen=1;
opts_knit.if_log=0;
%
results_knit=cell(1,nvariants);
warnings=[];
for ivariant=1:nvariants
    fullnames=cell(1,nsets);
    for iset=1:nsets
        fullnames{iset}=cat(2,coord_path,filesep,results_res{iset,ivariant}.coord_file);
    end
    aux=struct;
    aux.nsets=nsets;
    [data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
    disp('***********')
    disp(sprintf('variant %5.0f',ivariant));
    disp(sprintf('first file %s',fullnames{1}));
    disp(sprintf(' last file %s',fullnames{nsets}));
    %
    [data_aligned,aux_align_out]=rs_align_coordsets(data_read,setfields(struct(),{'opts_align','opts_check'},{opts_align,opts_check}));
    nstims_all=data_aligned.sets{1}.nstims;
    disp('aligned');
    %
    meth_string=results_res{iset,ivariant}.meth_string;
    imeth=find(contains(meth_strings,meth_string)>0);
    if length(imeth)==1
        opts_knit.allow_scale=if_allow_scales(imeth);
        aux.opts_knit=opts_knit;
        [data_consensus,aux_knit_out]=rs_knit_coordsets(data_aligned,setfields(struct(),{'opts_knit','opts_check'},{opts_knit,opts_check}));
        disp(sprintf('knit with dim_max_in=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',opts_knit.dim_max_in,opts_knit.pcon_init_method,opts_knit.allow_scale));
        if if_remove_pipeline
            data_consensus.sets{1}=struct();
        end
        results_knit{ivariant}.data_consensus=data_consensus;
        %collect key variance and shuffle stats
        results_knit{ivariant}.rmsavail_overall=aux_knit_out.knit_stats.rmsavail_overall;
        results_knit{ivariant}.rmsdev_overall=aux_knit_out.knit_stats.rmsdev_overall;
        results_knit{ivariant}.rmsdev_overall_shuff=reshape(aux_knit_out.knit_stats.rmsdev_overall_shuff(:,1,1,:,1),[dim_max_in nshuffs]); %d1: dimension, d2: which shuffle
    else
        wmsg=sprintf('for variant %4.0f, method %s not recognized',ivariant,meth_string);
        warning(wmsg);
        warnings=strvcat(warnings,wmsg);
    end
end
if ~isempty(warnings)
    disp(warnings);
end
%reorganize according to options for coordinate calculation
%d1: spatial filter, d2: df/f or z, d3: p_crit for pca filtering, d4: dim red method, d5: mean subtract or not
results_knit=reshape(results_knit,results_dims(2:end));
%
%plot
%
resp_measures={'deltaF/F','z'};
n_sfs=size(results_knit,1);
n_resps=size(results_knit,2);
n_pcrits=size(results_knit,3);
n_meths=size(results_knit,4);
n_sm=size(results_knit,5);
%
p_shuffle=getinp('p-value for showing significance of frac var explained by shuffle test','f',[0 1],0.05);
for rm_ptr=1:n_resps
    rm_string=resp_measures{rm_ptr};
    for submean=0:n_sm-1
        if submean
            sm_string='-sm';
        else
            sm_string='';
        end
        tstring=cat(2,rm_string,sm_string,sprintf(' nsets: %1.0f',nsets));
        figure;
        set(gcf,'Position',[50 50 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        for isf=1:n_sfs
            for pcrit_ptr=1:n_pcrits
                isub=isf+(pcrit_ptr-1)*n_sfs;
                subplot(n_pcrits,n_sfs,isub)
                rk=squeeze(results_knit(isf,rm_ptr,pcrit_ptr,:,1+submean)); %rk has all the dimension reduction methods
                %plot frac var explained, scale of [0 1]
                hp_leg=[];
                for k=1:n_meths
                    fvex=1-rk{k}.rmsdev_overall./rk{k}.rmsavail_overall;
                    pts_sel=[1:dim_max_in];
                    hp=plot(pts_sel,fvex(pts_sel),'k');
                    set(hp,'LineStyle','none'); %to add line if significant
                    set(hp,'Color',meth_colors{k});
                    set(hp,'Marker',meth_markers{k});
                    hold on;
                    hp_leg(k)=hp;
                    %what fraction of shuffles do better?
                    nshuffs_have=size(rk{k}.rmsdev_overall_shuff,2);
                    shuff_frac=sum(repmat(rk{k}.rmsdev_overall,1,nshuffs_have)>rk{k}.rmsdev_overall_shuff,2)/nshuffs_have;
                    %highlight the values with shuff_frac<p_shuffle
                    pts_sig=find(shuff_frac<=p_shuffle);
                    for p=1:length(pts_sig)
                        if pts_sig(p)>1
                            pts_sel=[-1 0]+pts_sig(p);
                            hs=plot(pts_sel,fvex(pts_sel),'k');
                            set(hs,'Color',meth_colors{k});
                            set(hs,'Marker',meth_markers{k});
                        end
                    end
                end
                if (isub==1)
                    legend(hp_leg,meth_legs,'FontSize',7,'Interpreter','none','Location','SouthWest');
                end
                set(gca,'XLim',[1 dim_max_in]);
                set(gca,'XTick',[1:dim_max_in]);
                set(gca,'YLim',[0 1]);
                xlabel('dim');
                ylabel('frac var expl');
                  %legends are meth_legs{:}
                title(sprintf('sf %1.0f pcrit %5.3f',sfilts(isf),pcrits(pcrit_ptr)));
            end
        end
        axes('Position',[0.01,0.01,0.01,0.01]);
        text(0,0,cat(2,tstring,' ',results_file),'Interpreter','none');
        axis off;
        axes('Position',[0.01,0.05,0.01,0.01]);
        text(0,0,sprintf('p_shuffle: %5.3f',p_shuffle),'Interpreter','none');
        axis off
        %
    end %submean
end %rm_ptr
clear results results_res aux_align_out aux_knit_out data_aligfned data_read rk
disp('suggest saving the workspace for re-use');
