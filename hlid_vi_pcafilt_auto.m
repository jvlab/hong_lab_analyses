%hlid_vi_pcafilt_auto: KC volumetric imaging, data from George Barnum, Hong Lab,
%
% derived from hlid_vi_pca_filt, automated, and can run over several datasets
% can customize pcrits, heatmap range (d_range_*)
%
% * no spatial heatmap, timecourse or Pearson plots
% * key variables saved in results
% * pca selection based on p-value of f-ratio of stim x repeat array
% 
%   See also:  HLID_VI_READ, HLID_VI_PCAFILT, HLID_VI_PCAFILT_SUMM, HLID_VI_SPATIALFILTER, HLID_VI_STIMNAMES, HLID_VARRATS, HLID_VI_EXPLORE, HLID_VI_PCASELECT.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_files') data_files={...
    '20240502_a_30s_output_walk_mc_2.hdf5',...
    '20240502_a_med_30s_output_walk_mc_2.hdf5',...
    '20241010_c_30s_output_walk.hdf5',...
    '20241010_c_30s_output_walk_mc2.hdf5',...
    '20241010_c_med_30s_output_mc_2.hdf5',...
    '20241015_a_30s_output_walk.hdf5',...
    '20241103_a_30s_output_walk.hdf5',...
    '20241230_a_med_30s_output_mc_2.hdf5',...
    '20241230_b_med_30s_output_mc_2.hdf5',...
    '20241027_a_30s_output_walk.hdf5',...
    '20241027_a_med_30s_output_mc_2.hdf5',...
    '20250807_b_2_30s_output_fixed.hdf5',...
    '20250808_b_30s_output_fixed.hdf5'};
end
%
if ~exist('d_range_avg') d_range_avg=[-1 1]; end %colormap range for average responses
if ~exist('d_range_indiv') d_range_indiv=d_range_avg; end %colormap range for average responses
%
if ~exist('n_stims') n_stims=24; end
if ~exist('n_repts') n_repts=5; end
if ~exist('max_timepoints') max_timepoints=0; end %until end of file
%
if ~exist('logrange') logrange=10^2; end %range of pc powers to plot
if ~exist('pcrits') pcrits=[1.0 0.2 0.1 0.05 0.025]; end %critical values of F to plot
%
dist_names={'cosine'};
n_dists=1;
%dist_names={'cosine','Pearson'};
%n_dists=length(dist_names);
%
resp_measures={'deltaF/F','z'};
for k=1:length(resp_measures)
    disp(sprintf('%1.0f->response measure %s',k,resp_measures{k}));
end
resp_measure=resp_measures{getinp('choice','d',[1 length(resp_measures)],1)};
%
if ~exist('opts_read') opts_read=struct;end
%
opts_read.max_timepoints=0; %read all time points
opts_read.if_remnan=1; %remove NaN's looking at all of the data
opts_read.if_spatialfilter=getinp('1 for spatial filter','d',[0 1],0);
if opts_read.if_spatialfilter
    opts_read.sfilt_hw=0.5*getinp('spatial kernel full width (0 is no filter)','d',[0 Inf],1);
    sf_string=sprintf('sf: kernel hw=%3.1f',opts_read.sfilt_hw);
else
    sf_string='sf: none';
end
%
if ~exist('if_timemark1') if_timemark1=0; end
if ~exist('if_timemark2') if_timemark2=1; end
%
if_spatial_pattern=0; %getinp('1 to plot spatial patterns','d',[0 1],1);
if_temporal_pattern=0; %getinp('1 to plot temporal patterns (-1: with extended baseline)','d',[-1 1],-1);
if_dist=1; % getinp('1 to plot distances','d',[0 1],1);
%
if_reorder=1; % getinp('1 to reorder stimuli','d',[0 1],1);
rept_list=[1:n_repts]; % getinp('repeat list','d',[1 n_repts],1:n_repts);
stim_list=[1:n_stims]; % getinp('stimulus list','d',[1 n_stims],1:n_stims);
%
opts_read.data_path=data_path;
%opts_read.data_file=data_file;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
if if_temporal_pattern==-1
    opts_read.if_keep_all_raw=1;
else
    opts_read.if_keep_all_raw=0;
end
%
n_pcasels=length(pcrits); %number of ways to select PC's
n_cols=n_pcasels+1; %first col for pc plots, rest for correlations
%
results=cell(1,n_files);
for file_ptr=1:n_files
    ifile=data_files_selected(file_ptr);
    data_file=data_files{ifile};
    results{file_ptr}.data_file=data_file;
    %
    disp('***********');
    disp(sprintf('processing file %3.0f of %3.0f: %s',file_ptr,n_files,data_file));
    opts_read.data_file=data_file;
    [s,opts_read_used]=hlid_vi_read(opts_read);
    n_pixels=s.n_pixels_kept;
    %
    read_data_file_short=strrep(s.opts_read.data_file,'.hdf5','');
    rept_string=cat(2,'repts: ',sprintf(' %2.0f',opts_read.rept_list));
    sf_tp_string=cat(2,sf_string,sprintf('; timepoints: [0 %2.0f]',s.n_timepoints_read));
    %
    tstring=cat(2,read_data_file_short,':',sf_tp_string,',',resp_measure,', ',rept_string);
    disp(sprintf('read %s',tstring));
    disp(s)
    %
    stims=hlid_vi_stimnames;
    %
    %compute mean response measure (delta-F/F or z), averaged over repeats
    %
    resp_maxlength=size(s.responses,2);
    deltaF=s.responses-repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
    switch resp_measure
        case 'deltaF/F'
            v=deltaF./repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
        case 'z'
            v=deltaF./repmat(reshape(s.baseline_stdvs,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 resp_maxlength 1 1]);
    end
    clear deltaF
    resp_minlength=sum(0==any(any(any(isnan(v),1),3),4));
    v_indiv_repts=reshape(v(:,[1:resp_minlength],:),[n_pixels*resp_minlength,s.n_repts_kept*s.n_stims_kept]);
    %
    if if_temporal_pattern==-1
        deltaF_alltime=s.pixel_values_kept-repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 s.n_timepoints_read 1 1]);
        switch resp_measure
            case 'deltaF/F'
                v_alltime=deltaF_alltime./repmat(reshape(s.baseline_means,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 s.n_timepoints_read 1 1]);
            case 'z'
                v_alltime=deltaF_alltime./repmat(reshape(s.baseline_stdvs,[s.n_pixels_kept,1,s.n_repts_kept,s.n_stims_kept]),[1 s.n_timepoints_read 1 1]);
        end
        clear deltaF_alltime
    end
    xyz_range=double([min(s.xyz_kept);max(s.xyz_kept)]);
    x_len=xyz_range(2,1)-xyz_range(1,1)+1;
    y_len=xyz_range(2,2)-xyz_range(1,2)+1;
    xyz_rel=double(s.xyz_kept)-repmat(xyz_range(1,:),s.n_pixels_kept,1)+1; %relative index of all kept pixels
    %
    display_ptr_order=[1:s.n_stims_kept];
    if if_reorder
        display_sort=zeros(1,s.n_stims_kept);
        for k=1:s.n_stims_kept
            display_sort(k)=find(stims.display_order==s.opts_read.stim_list(k));
        end
        [dsort,display_ptr_order]=sort(display_sort);
    end
    %
    %do pca and show properties
    %
    n_pc_max=s.n_repts_kept*s.n_stims_kept;
    [svd_u,svd_s,svd_v]=svd(v_indiv_repts,'econ'); %data=u*s*v'
    %
    disp('computed pcs')
    svd_v_max=max(abs(svd_v(:)));
    svd_u_max=max(abs(svd_u(:)));
    svd_s_dsq=diag(svd_s).^2;
    %
    frats=ones(n_pc_max,1);
    for ipc=1:n_pc_max
        spatem=reshape(svd_u(:,ipc),[n_pixels resp_minlength]);
        repstm=reshape(svd_v(:,ipc),[s.n_repts_kept s.n_stims_kept]);
        varrats=hlid_varrats(reshape(repstm,[1 s.n_repts_kept,s.n_stims_kept]));
        frats(ipc)=varrats.frat;
    end
    frats_pvals=1-fcdf(frats,varrats.fdof(1),varrats.fdof(2));
    %  
    hf=figure;
    set(gcf,'NumberTitle','off');
    set(gcf,'Position',[50 50 1400 800]);
    set(gcf,'Name',cat(2,' pca sel: ',tstring));
    %
    subplot(3,n_cols,1)
    semilogy(svd_s_dsq,'.-');
    xlabel('eigenvalue');
    set(gca,'XLim',[0 n_pc_max]);
    ylabel(cat(2,resp_measure,' var explained'));
    set(gca,'YLim',max(svd_s_dsq)*[1/logrange 1]);
    %
    subplot(3,n_cols,1+n_cols);
    semilogy(frats);
    hold on;
    for icrit=1:length(pcrits)
        finv_crit=finv(1-pcrits(icrit),varrats.fdof(1),varrats.fdof(2));
        semilogy([0 n_pc_max],repmat(finv_crit,1,2),'k:');
    end
    xlabel('eigenvalue');
    set(gca,'XLim',[0 n_pc_max]);
    ylabel('F ratio');
    set(gca,'YLim',[.1 100]);
    %
    subplot(3,n_cols,1+2*n_cols);
    semilogy(frats_pvals);
    hold on;
    for icrit=1:length(pcrits)
        semilogy([0 n_pc_max],repmat(pcrits(icrit),1,2),'k:');
    end
    xlabel('eigenvalue');
    set(gca,'XLim',[0 n_pc_max]);
    ylabel('p(F ratio)');
    set(gca,'YLim',[10^-5 1]);
    %
    %reconstruct with pcas selected based on F-ratio for stims x repeats
    %
    stim_names=stims.names_short(s.opts_read.stim_list(display_ptr_order));
    %
    results{file_ptr}.stim_names=stim_names;
    results{file_ptr}.dist_names=dist_names;
    results{file_ptr}.resp_measure=resp_measure;
    results{file_ptr}.eiv_squared=svd_s_dsq;
    %statistics for each pc
    results{file_ptr}.frats=frats;
    results{file_ptr}.frats_pvals=frats_pvals;
    results{file_ptr}.fdof=varrats.fdof;
    %
    for i_pcasel=1:n_pcasels
        %logic of hlid_vi_pcaselect, based on max p-value
        opts_pcasel=struct();
        opts_pcasel.n_pc_max=n_pc_max;
        opts_pcasel.eiv_squared=svd_s_dsq;
        opts_pcasel.frats=frats;
        opts_pcasel.fdof=varrats.fdof;
        pval_max=pcrits(i_pcasel);
        frat_min=finv(1-pval_max,opts_pcasel.fdof(1),opts_pcasel.fdof(2));
        pc_sel=find(opts_pcasel.frats>=frat_min);
        %
        frac_var_used=sum(opts_pcasel.eiv_squared(pc_sel))/sum(opts_pcasel.eiv_squared(:));
        pc_sel_string=sprintf('p %5.3f F %5.3f npc %3.0f pfrac %5.3f',pval_max,frat_min,length(pc_sel),frac_var_used);
        disp(sprintf('pca set %2.0f: %s',i_pcasel,pc_sel_string));
        %
        proj_pc=svd_v(:,pc_sel)*svd_v(:,pc_sel)'; %projection on selected pcs
        vfilt_indiv_repts=v_indiv_repts*proj_pc;
        vfilt=reshape(reshape(v,[n_pixels*resp_maxlength,s.n_repts_kept*s.n_stims_kept])*proj_pc,size(v));
        %
        results{file_ptr}.pcrit(i_pcasel)=pval_max;
        results{file_ptr}.frac_var_used(i_pcasel)=frac_var_used;
        results{file_ptr}.npc(i_pcasel)=length(pc_sel);
        %
        % average across repeats
        %
        v_across_repts=reshape(mean(vfilt,3,'omitnan'),[s.n_pixels_kept,resp_maxlength,s.n_stims_kept]);
        v_across_repts=reshape(v_across_repts(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_stims_kept]);
        %
        for idist=1:n_dists
            switch dist_names{idist}
                case 'cosine'
                    dot_prods=v_across_repts'*v_across_repts;
                    mags=sqrt(diag(dot_prods));
                    heatmap=dot_prods./(mags*mags');
                case 'Pearson'
                    heatmap=corr(v_across_repts);
            end
            subplot(2,n_cols,1+i_pcasel);
            heatmap_ordered=heatmap(display_ptr_order,display_ptr_order);
            results{file_ptr}.corrs_avg(:,:,i_pcasel,idist)=heatmap_ordered;
            imagesc(heatmap_ordered-diag(diag(heatmap_ordered)),d_range_avg); %remove the diagonal so as not to inflate the scale
            set(gca,'XTick',[1:s.n_stims_kept]);
            set(gca,'XTickLabels',stim_names);
            set(gca,'YTick',[1:s.n_stims_kept]);
            set(gca,'YTickLabels',stim_names);
            axis square;
            title(pc_sel_string,'FontSize',8,'FontWeight','normal');
        end
        clear v_across_repts
        %
        % keep individual repeats
        %
        v_indiv_repts_dist=reshape(vfilt,[s.n_pixels_kept,resp_maxlength,s.n_repts_kept*s.n_stims_kept]);
        v_indiv_repts_dist=reshape(v_indiv_repts_dist(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_repts_kept*s.n_stims_kept]);
        %
        var_rats=hlid_varrats(reshape(v_indiv_repts_dist,[s.n_pixels_kept*resp_minlength,s.n_repts_kept,s.n_stims_kept]));
        disp(var_rats);
        results{file_ptr}.var_rats{i_pcasel}=rmfield(var_rats,'centroid');
        %
        dpo_expanded=s.n_repts_kept*repmat(display_ptr_order-1,s.n_repts_kept,1); %expand the display pointer order to take into account blocks of repeats
        dpo_expanded=dpo_expanded+repmat([1:s.n_repts_kept]',1,s.n_stims_kept);
        dpo_expanded=dpo_expanded(:)';
        tick_locs=[1:s.n_stims_kept]*s.n_repts_kept-(s.n_repts_kept-1)/2; %tick locations for a block of repeats
        %
        for idist=1:n_dists
            switch dist_names{idist}
                case 'cosine'
                    dot_prods=v_indiv_repts_dist'*v_indiv_repts_dist;
                    mags=sqrt(diag(dot_prods));
                    heatmap=dot_prods./(mags*mags');
                case 'Pearson'
                    heatmap=corr(v_indiv_repts_dist);
            end
            subplot(2,n_cols,1+i_pcasel+(1+n_pcasels));
            heatmap_expanded=heatmap(dpo_expanded,dpo_expanded);
            results{file_ptr}.corrs_indiv(:,:,i_pcasel,idist)=heatmap_expanded;
            imagesc(heatmap_expanded-diag(diag(heatmap_expanded)),d_range_indiv); %remove the diagonal so as not to inflate the scale
            set(gca,'XTick',tick_locs);
            set(gca,'XTickLabels',stim_names);
            set(gca,'YTick',tick_locs);
            set(gca,'YTickLabels',stim_names);
            axis square;
            title(sprintf('overall F ratio: %8.4f',var_rats.frat),'FontSize',8,'FontWeight','normal');
        end
        %
    end %i_pcasel
    %
    axes('Position',[0.01,0.01,0.01,0.01]);
    text(0,0,tstring,'Interpreter','none');
    axis off
    axes('Position',[0.75 0.51 0.01 0.01]);
    text(0,0,sprintf('range: [%4.2f,%4.2f]',d_range_avg));
    axis off;
    axes('Position',[0.75 0.01 0.01 0.01]);
    text(0,0,sprintf('range: [%4.2f,%4.2f]',d_range_indiv));
    axis off;
    %
    clear s svd_u v v_indiv_repts v_indiv_repts_dist vfilt*
    %
end %file_ptr
