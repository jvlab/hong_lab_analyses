%hlid_rastim_mds_coords_make: read calcium imaging data from Hong Lab, 
% create coordinate sets from several transformations of the z-scores (averaged across datasets)
% and create consensus with jackknifes on stimuli and datasets
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
% differences c/w hlid_rastim_mds_coords_demo:
%   * does basic variance statistics on consensus, but not shuffle-based statistics
%   * no plotting
%   * number of coords limited to number of stimuli
%   * pipeline only saved for non-jackknifed data and non-jackknifed knit
%   * max dimension across methods asserted to be nstims-1, reducing by 1 because of jackknife on stimuli
%   * all datasets must have all stimuli (to allow for jackknifing on stimuli)
%   * uses stim_labels [short names of stimuli (omitting concentration)] for output files
%   * eigenvalue statistics for individual files not saved when jackknifing on files (since these are computed within files)
%
% results saved in r_all, with subfields for all data, and the two kinds of jackknifes
%
%  See also:  HLID_LOCALOPTS, HLID_READ_COORDDATA_DEMO, HLID_DA_STIMSELECT, HLID_RASTIM2COORDS, DOMDS, PSG_ISOMAP_DEMO
%   PSG_KNIT_STATS, HLID_RASTIM_MDS_COORDS_SUMM, HLID_RASTIM_MDS_COORDS_SUMM2, LID_RASTIM_MDS_COORDS_SUBSAMP,
%   HLID_RASTIM_MDS_COORDS_DEMO.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
nshuffs=0; % no statistics
if_debug=getinp('1 for debugging','d',[0 1]);
if_keep_jackfile_data=getinp('1 to keep data fields when jackknifing on files (0 to save space)','d',[0 1],0);
%
if ~exist('neigs_show_max') neigs_show_max=10; end
if ~exist('dim_max_in_def') dim_max_in_def=10; end %maximum dimension for calculating consensus
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
% analyze with and without subtracting the mean
%
nsubs=1+if_submean; 
maxdim_all=min([min(nrois),min(nstims_each)])-if_submean;
%
dim_max_in=getinp('max dim for consensus statistics','d',[1 nstims-1],dim_max_in_def); %nstims-1 since we will jackknife on stimuli
%
r_all=struct;
for jack_mode=1:3 %all data, and jackknife on files, and jackknife on stimuli
    switch jack_mode
        case 1 %all files, all stimuli
            jack=0;
            njack=1;
            r_field='jackknife_omit_none';
        case 2 %jackknife on files, all stimuli
            jack=nfiles;
            njack=nfiles;
            r_field='jackknife_by_file';
        case 3 % all files, jackknife on stimuli
            jack=nstims;
            njack=nstims;
            r_field='jackknife_by_stim';
    end
    if (if_debug)
        njack=min(njack,2);
    end
    r_all.(r_field)=cell(1,njack);
    for ijack=1:njack
        disp('*********');
        disp(sprintf('processing jack_mode %1.0f jackknife %2.0f of %2.0f',jack_mode,ijack,njack));
        r=struct;
        if jack_mode~=2 %these are computed within datasets, so no need to recompute if jackknifed on datasets
            r.eivals=cell(nfiles,nmeths,nsubs); %files, methods, sub mean?
            r.coords=cell(nfiles,nmeths,nsubs);
            r.maxdim=zeros(nfiles,nsubs);
            r.pwr_ratio=zeros(nfiles,nmeths,nsubs);
            r.part_ratio=zeros(nfiles,nmeths,nsubs);
        end
        %
        max_dev=zeros(nfiles,nsubs);
        svd_dev=zeros(nfiles);
        %implement jackknifing
        stims_keep=[1:nstims];
        files_keep=[1:nfiles];
        file_removed=0;
        stim_removed=0;
        switch jack_mode
            case 1
                jack_mode_desc='jackknife: none';
            case 2
                jack_mode_desc='jackknife: on files';
                files_keep=setdiff([1:nfiles],ijack);
                file_removed=ijack;
            case 3
                jack_mode_desc='jackknife: on stimuli';
                stims_keep=setdiff([1:nstims],ijack);
                stim_removed=istim;
        end
        %
        nfiles_keep=length(files_keep);
        nstims_keep=length(stims_keep);
        for ifile_ptr=1:nfiles_keep
            ifile=files_keep(ifile_ptr);
            for submean=0:if_submean
                ru=resps_use{ifile};
                ru=ru(stims_keep,:);
                if (submean)
                    ru=ru-repmat(mean(ru,1),nstims_keep,1);
                end
                r.maxdim(ifile,1+submean)=min(nrois(ifile),nstims_keep)-submean;
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
                   if size(coords,2)>size(coords,1)
                       coords=coords(:,1:size(coords,1));
                   end
                   if jack_mode~=2
                        r.eivals{ifile,imeth,1+submean}=eivals;
                        r.pwr_ratio(ifile,imeth,1+submean)=sum(eivals(eivals>0))/sum(abs(eivals));
                        r.part_ratio(ifile,imeth,1+submean)=(sum(abs(eivals)).^2)/sum(eivals.^2);
                   end
                    r.coords{ifile,imeth,1+submean}=coords;
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
        end %ifile_ptr
        %save file names, convert to das, sas, sets structures
        r.nstims=nstims_keep;
        r.nrois=nrois(files_keep);
        r.meths=meths;
        r.filenames_short=filenames_short(files_keep);
        r.data_desc='d1: method, d2: mean not subtracted, then mean subtracted';
        r.data=cell(nmeths,1+if_submean);
        for imeth=1:nmeths
            for submean=0:if_submean
                %use nfiles_keep and ifile_ptr for r.data, so that there are no missing entries
                r.data{imeth,1+submean}.ds=cell(nfiles_keep,1);
                r.data{imeth,1+submean}.sas=cell(nfiles_keep,1);
                r.data{imeth,1+submean}.sets=cell(nfiles_keep,1);
                for ifile_ptr=1:nfiles_keep
                    ifile=files_keep(ifile_ptr);
                    coords=r.coords{ifile,imeth,1+submean};
                    nd=min(nstims_keep,size(coords,2));
                    r.data{imeth,1+submean}.ds{ifile_ptr}=cell(1,nd);
                    for k=1:nd
                        r.data{imeth,1+submean}.ds{ifile_ptr}{k}=coords(:,1:k);
                    end %dim
                    r.data{imeth,1+submean}.sas{ifile_ptr}=struct;
                    r.data{imeth,1+submean}.sas{ifile_ptr}.nstims=nstims_keep;
 %                   typenames_all=das{ifile}.response_amplitude_stim.stim(stims_nonan{ifile})';
                    r.data{imeth,1+submean}.sas{ifile_ptr}.typenames=stim_labels(stims_keep); %use short names of stimuli
                    btc_specoords_all=eye(nstims);
                    r.data{imeth,1+submean}.sas{ifile_ptr}.btc_specoords=btc_specoords_all(stims_keep,:);
                    %only keep the stimuli that are kept by jackknife
                    r.data{imeth,1+submean}.sets{ifile_ptr}=struct;
                    r.data{imeth,1+submean}.sets{ifile_ptr}.type='data';
                    r.data{imeth,1+submean}.sets{ifile_ptr}.dim_list=[1:nd];
                    r.data{imeth,1+submean}.sets{ifile_ptr}.nstims=nstims_keep;
                    r.data{imeth,1+submean}.sets{ifile_ptr}.label_long=cat(2,pathname,filesep,filenames_short{ifile});
                    r.data{imeth,1+submean}.sets{ifile_ptr}.label=filenames_short{ifile};
                    if (jack_mode==1)
                        r.data{imeth,1+submean}.sets{ifile_ptr}.pipeline.type='mds';
                        r.data{imeth,1+submean}.sets{ifile_ptr}.pipeline.opts=meths{imeth};
                        r.data{imeth,1+submean}.sets{ifile_ptr}.pipeline.opts.if_submean=submean;
                    else
                        r.data{imeth,1+submean}.sets{ifile_ptr}.pipeline=[];
                    end
                end %ifile_ptr
              end %submean
        end %imeth
        %
        %align and knit with rs package
        aux=struct;
        aux.opts_align.if_log=0;
        aux.opts_knit.if_log=0;
        aux.opts_knit.allow_scale=0;
        aux.opts_knit.if_normscale=0;
        aux.opts_knit.keep_details=0; 
        aux.opts_knit.if_stats=1; %do basic stats but no shuffles
        aux.opts_knit.nshuffs=0;
        aux.opts_knit.if_plot=0; %no plots
        aux.opts_knit.dim_max_in=dim_max_in;
        %
        r.data_knit=cell(nmeths,1+if_submean);
        r.knit_stats=cell(nmeths,1+if_submean);
        %
        disp(' ');
        for imeth=1:nmeths
            for submean=0:if_submean
                meth_text=cat(2,meths{imeth}.name_short,sprintf(' submean=%1.0f',submean));
                disp(sprintf('processing %s, jack_mode %1.0f jack %2.0f',meth_text,jack_mode,ijack));
                data_in=r.data{imeth,1+submean};
                [data_align,aux_align]=rs_align_coordsets(data_in,aux);
                [data_knit,aux_knit]=rs_knit_coordsets(data_align,aux);
                if (jack_mode>1)
                    data_knit=rmfield(data_knit,'sets');
                    data_knit=rmfield(data_knit,'sas');
                end
                r.data_knit{imeth,1+submean}=data_knit;
                r.knit_stats{imeth,1+submean}=aux_knit.knit_stats;

            end %submean
        end %imeth
        %document the jackknifing
        r.jack_mode_desc=jack_mode_desc;
        r.jack_mode=jack_mode;
        r.jack_count=njack;
        if jack_mode==2
            r.file_removed=filenames_short{ijack};
        else
            r.file_removed=0;
        end
        if jack_mode==3
            r.stim_removed=stim_labels{ijack};
        else
            r.stim_removed=0;
        end
        if ((if_keep_jackfile_data==0) & (jack_mode==2));
            r=rmfield(r,'data');
            r=rmfield(r,'coords');
        end
        r_all.(r_field){ijack}=r; %save jackknifed values
    end %ijack
end %jack_mode (standard or jackknife)
