%hlid_rastim_mds_coords_demo: read calcium imaging data from Hong Lab, 
% examine several transformations of the z-scores
% 
% variants include:
%   Euclidean distances between z-scored values, with dimension reduction
%     via svd of z-scored values, or mds of distances
%   angular distance (arccos) between dot products of normalized z-scored values 
%   arc length distance between normalized z-scored values
%   angular distance (arccos) between dot products of normalized z-scored values after subtraction of mean across ROIs (Pearson)
%   arc length distance between normalized z-scored values ater subtraction of mean across ROIs (Pearson)
%
% Number of ROIs without NaNs must be at least as large as number of stimuli.
%
% Note that centering (subtraction of mean across ROI's) is not the same as "subtract the mean"
%  option in hlid_rastim2coords, which subtracts the mean across stimuli
%
% results saved in r
%
% 25Oct25: begin saving sa, ds, sets so that a modularized version of psg_knit_stats_demo, psg_align_vara_demo will run
%
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, HLID_DA_STIMSELECT, HLID_RASTIM2COORDS, DOMDS, PSG_ISOMAP_DEMO
%   PSG_KNIT_STATS, HLID_RASTIM_MDS_COORDS_SUMM.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
if ~exist('neigs_show_max') neigs_show_max=10; end
%
if ~exist('opts_dasel')
    opts_dasel=struct;
end
%
if ~exist('meths')
    meths=cell(0);
    meths{1}.name_full='Euclidean distance via SVD';
    meths{1}.name_short='Euc dist SVD';
    meths{1}.xform='none';
    meths{1}.dimred='svd';
    meths{1}.name_file='euc_svd';
    %
    meths{2}.name_full='Euclidean distance via MDS';
    meths{2}.name_short='Euc dist MDS';
    meths{2}.xform='none';
    meths{2}.name_file='euc_mds';
    %
    meths{3}.name_full='cosine similarity'; %1-normalized dot product
    meths{3}.name_short='cos sim';
    meths{3}.xform='1-dp';
    meths{3}.name_file='cos_sim';
    %
    meths{4}.name_full='cosine similarity as angle';
    meths{4}.name_short='cos ang';
    meths{4}.xform='acos(dp)';
    meths{4}.name_file='cos_ang';
    %
    meths{5}.name_full='cosine similarity as chord';
    meths{5}.name_short='cos chord';
    meths{5}.xform='sqrt(2)*sqrt(1-dp)';
    meths{5}.name_file='cos_chord';
    %
    meths{6}.name_full='Pearson similarity'; %1-normalized centered dot product
    meths{6}.name_short='Pearson sim';
    meths{6}.xform='1-dp';
    meths{6}.name_file='pears_sim';
    %
    meths{7}.name_full='Pearson similarity as angle';
    meths{7}.name_short='Pearson ang';
    meths{7}.xform='acos(dp)';
    meths{7}.name_file='pears_ang';
    %
    meths{8}.name_full='Pearson similarity as chord';
    meths{8}.name_short='Pearson chord';
    meths{8}.xform='sqrt(2)*sqrt(1-dp)';
    meths{8}.name_file='pears_chord';
%
end
nmeths=length(meths);
meths_allnames=cell(1,nmeths);
for imeth=1:nmeths
    if ~isfield(meths{imeth},'dimred')
        meths{imeth}.dimred='mds';
    end
    meths_allnames{imeth}=meths{imeth}.name_short;
end
%
if ~exist('HongLab_fn') HongLab_fn='C:/Users/jdvicto/Dropbox/From_HongLab/HongLabOrig_for_jdv/data/kc_soma_nls/2022-10-10__fly01__megamat0.mat'; end
%
ui_prompt='Select one or more raw data files';
ui_filter='*.mat';
if_ok=0;
while (if_ok==0)
    [filenames_short,pathname,filter_index]=uigetfile(ui_filter,ui_prompt,'Multiselect','on');
    if ~iscell(filenames_short)
        filenames_short={filenames_short};
    end
    nfiles=length(filenames_short);
    % for ifile=1:nfiles
    %     disp(sprintf(' file %1.0f->%s',ifile,filenames_short{ifile}));
    % end
    das=cell(1,nfiles);
    opts_used_daset=cell(1,nfiles);
    dsids=cell(1,nfiles);
    if_ok=1;
    for ifile=1:nfiles
        %
        da=load(cat(2,pathname,filesep,filenames_short{ifile}));
        dsid=da.meta.title;
        %'-'can be used within fields of file name
        dsid=strrep(dsid,'/','-');
        dsid=strrep(dsid,'\','-');
        dsid=strrep(dsid,'_','-');
        while contains(dsid,'--')
            dsid=strrep(dsid,'--','-');
        end
        if (ifile==1)
            stimulus_names=strvcat(da.response_amplitude_stim.stim');
            nstims=length(stimulus_names);
        end
        stimulus_names_this=strvcat(da.response_amplitude_stim.stim');
        nstims_this=length(stimulus_names_this);
        if_match=1;
        if nstims~=nstims_this
            if_match=0;
        else
            if any(stimulus_names~=stimulus_names_this)
                if_match=0;
            end
        end
        dsids{ifile}=dsid;
        %
        [das{ifile},optsused_dasel{ifile}]=hlid_da_stimselect(da,opts_dasel);
        disp(sprintf(' file %1.0f->%30s->%s',ifile,filenames_short{ifile},dsids{ifile}));
        if (if_match==0)
            disp(sprintf('file %2.0f: incompatible number of stimuli or stimulus names',ifile));
            if_ok=0;
        end
    end
    disp(sprintf(' %3.0f stimuli',nstims));
    if (if_ok==1)
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
%
resps_use=cell(1,nfiles);
stims_nonan=cell(1,nfiles);
nrois=zeros(1,nfiles);
nstims_each=zeros(1,nfiles);
for ifile=1:nfiles
    resps_raw=das{ifile}.response_amplitude_stim.mean_peak;
    nancols=find(all(isnan(resps_raw),1));
    resps=resps_raw(:,setdiff([1:size(resps_raw,2)],nancols));
    if size(resps,1)~=nstims
        warning(sprintf('number of responses (%2.0f) is not equal to number of stimuli (%2.0f)',size(resps,1),nstims));
    end
    if size(resps,1)>size(resps,2)
        warning(sprintf('number of responses (%2.0f) is greater than number of neural data channels or rois (%4.0f)',size(resps,1),size(resps,2)));
    end
    %
    stims_nonan{ifile}=find(all(~isnan(resps),2));
    resps_use{ifile}=resps(stims_nonan{ifile},:);
    nstims_each(ifile)=size(resps_use{ifile},1);
    nrois(ifile)=size(resps_use{ifile},2);
    %
    disp(sprintf('file %2.0f: ROIs found: %4.0f, ROIs that are all NaNs: %4.0f, ROIs analyzed: %4.0f,  number of stimuli without NaN responses: %3.0f',...
        ifile,size(resps_raw,2),length(nancols),nrois(ifile),nstims_each(ifile)));
end
%
stim_labels=das{1}.response_amplitude_stim.stim'; %stimulus names as a cell array
for istim=1:nstims
    if contains(stim_labels{istim},'@')
        stim_labels{istim}=stim_labels{istim}(1:min(find(stim_labels{istim}=='@')-1));
    end
    stim_labels{istim}=deblank(stim_labels{istim});
end
%if_submean=getinp('1 to also analyze with mean across stims subtracted','d',[0 1],1);
if_submean=1;
%
%
% analyze with and without subtracting the mean
%
nsubs=1+if_submean; 
maxdim_all=min([min(nrois),min(nstims_each)])-if_submean;
%
if_write=getinp('1 to write the consensus files','d',[0 1]);
if ~exist('write_prefix') write_prefix='hlid_consensus_coords'; end %must conrain hlid*_coords
if if_write
    filenames_out=cell(nmeths,nsubs);
    if_ok=0;
    while (if_ok==0)
        write_prefix=getinp('prefix for file output','s',[],write_prefix);
        disp(sprintf('data path is %s',pathname))
        write_suffix=getinp('suffix for file output, should designate dataset','s',[]);
        for imeth=1:nmeths
            for submean=0:if_submean
                filenames_out{imeth,1+submean}=cat(2,write_prefix,'_',write_suffix,'_',meths{imeth}.name_file);
                if submean
                    filenames_out{imeth,1+submean}=cat(2,filenames_out{imeth,1+submean},'-sm');
                end
                filenames_out{imeth,1+submean}=cat(2,filenames_out{imeth,1+submean},'.mat');
                disp(sprintf(' method %2.0f (%30s), submean=%1.0f: output name will be %s',imeth,meths{imeth}.name_full,submean,filenames_out{imeth,1+submean}));
            end
        end
        if_ok=getinp('1 if ok','d',[0 1]);
    end
end
%
r=struct;
r.eivals=cell(nfiles,nmeths,nsubs); %files, methods, sub mean?
r.coords=cell(nfiles,nmeths,nsubs);
r.maxdim=zeros(nfiles,nsubs);
r.pwr_ratio=zeros(nfiles,nmeths,nsubs);
r.part_ratio=zeros(nfiles,nmeths,nsubs);
%
max_dev=zeros(nfiles,nsubs);
svd_dev=zeros(nfiles);
for ifile=1:nfiles
    for submean=0:if_submean
        ru=resps_use{ifile};
        if (submean)
            ru=ru-repmat(mean(ru,1),nstims_each(ifile),1);
        end
        r.maxdim(ifile,1+submean)=min(nrois(ifile),nstims_each(ifile))-submean;
        ru_norm=ru./repmat(sqrt(sum(ru.^2,2)),1,nrois(ifile));
        ru_centered=ru-repmat(mean(ru,2),1,nrois(ifile));
        ru_centered_norm=ru_centered./repmat(sqrt(sum(ru_centered.^2,2)),1,nrois(ifile));
        [u_svd,s_svd,v_svd]=svd(ru); %for comparison with MDS of Euclidean distance
        eig_svd=sort(diag(s_svd),'descend');
        eig_svd=eig_svd(1:r.maxdim(ifile,1+submean));
        eig_svd=eig_svd/eig_svd(1);
        for imeth=1:nmeths
           %compute distance based on transformation method
           switch meths{imeth}.name_full
               case {'Euclidean distance via SVD','Euclidean distance via MDS'}
                   dists=sqrt(cootodsq(ru));
               case {'cosine similarity','cosine similarity as angle','cosine similarity as chord'}
                   dp=ru_norm*ru_norm';
               case {'Pearson similarity','Pearson similarity as angle','Pearson similarity as chord'}
                   dp=ru_centered_norm*ru_centered_norm';
               otherwise
                   disp(sprintf('distance %s not recognized',meths{imeth}.name_full));
           end
           switch meths{imeth}.xform
               case 'none'
               case '1-dp'
                   dists=1-dp;
               case 'acos(dp)'
                   dists=acos(dp);
               case 'sqrt(2)*sqrt(1-dp)'
                   dists=sqrt(2)*sqrt(1-dp);
            otherwise
                disp(sprintf('transform %s not recognized',meths{imeth}.xform));
           end
           %normalize all distances to rms of 1
           dnorm=sqrt(mean(dists(:).^2));
           dists=dists/dnorm;
           switch meths{imeth}.dimred
               case 'mds'
                   %code borrowed from psg_isomap
                    [eivals_raw,eivecs_raw]=domds(dists,1); % [eival,eivec]=domds(distmtx,p) does a multidimensional scaling
                    [eivals_sort,ix]=sort(abs(eivals_raw),'descend');
                    eivals_raw=eivals_raw(ix);
                    eivecs_raw=eivecs_raw(:,ix);
                    eivals=real(eivals_raw)/2; %per do_mds, the eigenvalues need to be divided by 2, prior to square root, to obtain coords that recapitulate the distances.
                    coords=eivecs_raw.*repmat(sqrt(abs(eivals')),[size(eivecs_raw,1),1]); %scale by eivals
               case 'svd'
                   [u,s,v]=svd(ru/dnorm);
                   sdiag_sq=(diag(s)).^2; %square to match mds convention
                   [eivals_sort,ix]=sort(abs(sdiag_sq),'descend');
                   eivals=sdiag_sq(ix);
                   coords=u*s;
           end
            r.eivals{ifile,imeth,1+submean}=eivals;
            r.coords{ifile,imeth,1+submean}=coords;
            r.pwr_ratio(ifile,imeth,1+submean)=sum(eivals(eivals>0))/sum(abs(eivals));
            r.part_ratio(ifile,imeth,1+submean)=(sum(abs(eivals)).^2)/sum(eivals.^2);
            %checks
            if all(eivals>0)
                dists_check=sqrt(cootodsq(coords));
                dev=max(abs(dists(:)-dists_check(:)));
                max_dev(ifile,1+submean)=max(max_dev(ifile,1+submean),dev);
            end
            if strcmp(meths{imeth}.name_full,'Euclidean distance via MDS')
                nc=min(length(eivals),length(eig_svd));
                svd_dev(ifile,1+submean)=max(abs(eivals(1:nc)/eivals(1)-eig_svd(1:nc).^2));
            end
        end %methods
        disp(sprintf(' for file %2.0f, submean %1.0f: max dev betw dists before and after MDS: %15.12f, max dev betw eivs from MDS and SVD (s.b. 0 if submean=1): %15.12f',...
            ifile,submean,max_dev(ifile,1+submean),svd_dev(ifile,1+submean)));
   end %subtract mean?
end %file
%
%plots
%
neigs_show=min(neigs_show_max,maxdim_all);
rng('default');
colors=rand(3,nfiles)'*0.75;
%
labels_short_nomm=cell(1,nfiles);
for ifile=1:nfiles
    last_sep=max([0,max(find(filenames_short{ifile}=='/')),max(find(filenames_short{ifile}=='\'))]);
    labels_short_nomm{ifile}=filenames_short{ifile}(last_sep+1:end);
    labels_short_nomm{ifile}=strrep(labels_short_nomm{ifile},'.mat','');
%    labels_short_nomm{ifile}=strrep(labels_short_nomm{ifile},'202','2');
%    labels_short_nomm{ifile}=strrep(labels_short_nomm{ifile},'fly','');
%    labels_short_nomm{ifile}=strrep(labels_short_nomm{ifile},'__','-');
end
labels_short=labels_short_nomm;
labels_short{ifile+1}='mean';
labels_short{ifile+2}='median';
%
for submean=0:if_submean
    figure;
    tstring=sprintf('submean=%1.0f',submean);
    set(gcf,'Position',[100 100 1400 800]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',tstring);
    %
    for imeth=1:nmeths
        subplot(2,nmeths,imeth);
        eiv_plot=NaN(neigs_show,nfiles);
        for ifile=1:nfiles
            neiv=min(neigs_show,r.maxdim(ifile,1+submean));
            eiv_plot(1:neiv,ifile)=r.eivals{ifile,imeth,1+submean}(1:neiv)./max(r.eivals{ifile,imeth,1+submean});
            hp=plot([1:neiv],eiv_plot(:,ifile),'k','LineWidth',1);
            set(hp,'Color',colors(ifile,:));
            hold on;
        end
        plot([1:neiv],mean(eiv_plot,2,'OmitNan'),'k','LineWidth',2);
        plot([1:neiv],median(eiv_plot,2,'OmitNan'),'k:','LineWidth',2);
        set(gca,'XLim',[1 neigs_show]);
        set(gca,'XTick',[1:neigs_show]);
        xlabel('eiv');
        set(gca,'YLim',[-.25 1]);
        title(meths{imeth}.name_short);
    end
    %
    %power ratio and participation ratio
    %
    for ipp=1:2
        switch ipp
            case 1
                vname='power ratio';
                v=r.pwr_ratio(:,:,1+submean);
            case 2
                vname='participation ratio';
                v=r.part_ratio(:,:,1+submean);
        end
        subplot(2,2,2+ipp);
        for ifile=1:nfiles
            hp=plot([1:nmeths],v(ifile,:),'k','LineWidth',1);
            set(hp,'Color',colors(ifile,:));
            hold on;
        end
        plot([1:nmeths],mean(v,1),'k','LineWidth',2);
        plot([1:nmeths],median(v,1),'k:','LineWidth',2);
        set(gca,'XTick',[1:nmeths]);
        set(gca,'XTickLabel',meths_allnames);
        if (ipp==1)
            set(gca,'YLim',[0.5 1]);
            legend(labels_short,'FontSize',6,'Interpreter','none','Location','SouthWest');
        else
            set(gca,'YLim',[0 max(get(gca,'YLim'))]);
        end
        title(vname);
    end
    axes('Position',[0.01,0.04,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none','FontSize',8);
    axis off;
end
%save file names, convert to das, sas, sets structures
r.nstims=nstims;
r.nstims_each=nstims_each;
r.nrois=nrois;
r.meths=meths;
r.filenames_short=filenames_short;
r.data_desc='d1: method, d2: mean not subtracted, then mean subtracted';
r.data=cell(nmeths,1+if_submean);
nd_max=Inf; %to find max dimension available, across all methods
for imeth=1:nmeths
    for submean=0:if_submean
        r.data{imeth,1+submean}.ds=cell(nfiles,1);
        r.data{imeth,1+submean}.sas=cell(nfiles,1);
        r.data{imeth,1+submean}.sets=cell(nfiles,1);
        for ifile=1:nfiles
            coords=r.coords{ifile,imeth,1+submean};
            nd=min(nstims_each(ifile),size(coords,2));
            nd_max=min(nd_max,nd);
            r.data{imeth,1+submean}.ds{ifile}=cell(1,nd);
            for k=1:nd
                r.data{imeth,1+submean}.ds{ifile}{k}=coords(:,1:k);
            end %dim
            r.data{imeth,1+submean}.sas{ifile}=struct;
            r.data{imeth,1+submean}.sas{ifile}.nstims=nstims_each(ifile);
            r.data{imeth,1+submean}.sas{ifile}.btc_specoords=eye(nd);
            %only keep the stimuli that are not NaN
            r.data{imeth,1+submean}.sas{ifile}.typenames=das{ifile}.response_amplitude_stim.stim(stims_nonan{ifile})';
            r.data{imeth,1+submean}.sets{ifile}=struct;
            r.data{imeth,1+submean}.sets{ifile}.type='data';
            r.data{imeth,1+submean}.sets{ifile}.dim_list=[1:nd];
            r.data{imeth,1+submean}.sets{ifile}.nstims=nstims_each(ifile);
            r.data{imeth,1+submean}.sets{ifile}.label_long=cat(2,pathname,filesep,filenames_short{ifile});
            r.data{imeth,1+submean}.sets{ifile}.label=filenames_short{ifile};
            r.data{imeth,1+submean}.sets{ifile}.pipeline.type='mds';
            r.data{imeth,1+submean}.sets{ifile}.pipeline.opts=meths{imeth};
            r.data{imeth,1+submean}.sets{ifile}.pipeline.opts.if_submean=submean;
        end %ifile
    end %submean
end %imeth
%
nshuffs=getinp('number of shuffles (0 for no stats)','d',[0 1000]);
dim_max_in=getinp('max dim for consensus statistics','d',[1 nd_max]);
%align and knit with rs package
aux=struct;
aux.opts_align.if_log=0;
aux.opts_knit.if_log=1;
aux.opts_knit.allow_scale=0;
aux.opts_knit.if_normscale=0;
aux.opts_knit.keep_details=0; 
aux.opts_knit.if_stats=double(nshuffs>0);
aux.opts_knit.nshuffs=nshuffs;
aux.opts_knit.if_plot=0; %plot it here
aux.opts_knit.dim_max_in=dim_max_in;
%
r.data_knit=cell(nmeths,1+if_submean);
r.knit_stats=cell(nmeths,1+if_submean);
%
for imeth=1:nmeths
    for submean=0:if_submean
        meth_text=cat(2,meths{imeth}.name_short,sprintf(' submean=%1.0f',submean))
        disp(' ');
        disp(sprintf('processing %s',meth_text));
        disp(' ');
        data_in=r.data{imeth,1+submean};
        [data_align,aux_align]=rs_align_coordsets(data_in,aux);
        [data_knit,aux_knit]=rs_knit_coordsets(data_align,aux);
        if (if_write)
            data_knit.sets{1}.pipeline.embedding_method=setfield(meths{imeth},'if_submean',submean);
            rs_write_coorddata(filenames_out{imeth,1+submean},data_knit);
        end
        r.data_knit{imeth,1+submean}=data_knit;
        %do stats if nshuffs>0
        if (nshuffs>0)
            knit_stats{imeth,1+submean}=aux_knit.knit_stats;
            knit_stats_setup=aux_knit.knit_stats_setup;
            knit_stats_setup.dataset_labels=strrep(filenames_short,'.mat','');
            knit_stats_setup.stimulus_labels=stim_labels;
            psg_knit_stats_plot(knit_stats{imeth,1+submean},knit_stats_setup);
            axes('Position',[0.5,0.07,0.01,0.01]); %for text
            text(0,0,meth_text,'Interpreter','none','FontSize',8);
            axis off;
            set(gcf,'Name',sprintf('consensus stats: %s',meth_text));
        end
    end %submean
end %imeth
