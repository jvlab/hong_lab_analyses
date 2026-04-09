%hlid_vi_pcafilt: KC volumetric imaging, data from George Barnum, Hong Lab,
%
%***to do* better eigenvalue selection
%
% derived from hlid_vi_explore, but allows filtering by characteristics of pcs
%
% * if_remnan forced to 1
% * max timepoints defaults to 93 (end of stimulus)
% * default set to read all data
% * defaults set to enable all plots
% * for extended time plots, all stimuli handled together, since this is needed for reconstruction from selected pcs
% 
%   See also:  HLID_VI_READ, HLID_VI_SPATIALFILTER, HLID_VI_STIMNAMES, HLID_VARRATS, HOLID_VI_EXPLORE.
%
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('data_file') data_file='20241027_a_30s_output_walk.hdf5'; end
if ~exist('n_stims') n_stims=24; end
if ~exist('n_repts') n_repts=5; end
if ~exist('max_timepoints') max_timepoints=0; end %until end of file
%
if ~exist('logrange') logrange=10^2; end %range of pc powers to plot
if ~exist('pcrits') pcrits=[0.001 0.01 0.05]; end %critical values of F to plot
%
dist_names={'cosine','Pearson'};
n_dists=length(dist_names);
%
resp_measures={'deltaF/F','z'};
for k=1:length(resp_measures)
    disp(sprintf('%1.0f->response measure %s',k,resp_measures{k}));
end
resp_measure=resp_measures{getinp('choice','d',[1 length(resp_measures)],1)};
%
if ~exist('opts_read') opts_read=struct;end
%
opts_read.max_timepoints=getinp('max timepoints to read (0: all)','d',[0 Inf],max_timepoints);
opts_read.if_remnan=1;
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
if_spatial_pattern=getinp('1 to plot spatial patterns','d',[0 1],1);
if_temporal_pattern=getinp('1 to plot temporal patterns (-1: with extended baseline)','d',[-1 1],-1);
if_dist=getinp('1 to plot distances','d',[0 1],1);
%
if_reorder=getinp('1 to reorder stimuli','d',[0 1],1);
rept_list=getinp('repeat list','d',[1 n_repts],1:n_repts);
stim_list=getinp('stimulus list','d',[1 n_stims],1:n_stims);
%
opts_read.data_path=data_path;
opts_read.data_file=data_file;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
if if_temporal_pattern==-1
    opts_read.if_keep_all_raw=1;
else
    opts_read.if_keep_all_raw=0;
end
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
%plot mean response (delta-F/F), averaged over repeats
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
[nrows,ncols]=nicesubp(s.n_stims_kept,0.7);
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
frats=ones(1,n_pc_max);
for ipc=1:n_pc_max
    spatem=reshape(svd_u(:,ipc),[n_pixels resp_minlength]);
    repstm=reshape(svd_v(:,ipc),[s.n_repts_kept s.n_stims_kept]);
    varrats=hlid_varrats(reshape(repstm,[1 s.n_repts_kept,s.n_stims_kept]));
    frats(ipc)=varrats.frat;
end
frats_pvals=1-fcdf(frats,varrats.fdof(1),varrats.fdof(2));
%  
figure;
set(gcf,'NumberTitle','off');
set(gcf,'Position',[50 50 1200 800]);
set(gcf,'Name',cat(2,sprintf(' eigenvalues: %2.0f, ',plane),tstring));
%
subplot(2,3,1)
semilogy(svd_s_dsq,'.-');
xlabel('eigenvalue');
set(gca,'XLim',[0 n_pc_max]);
ylabel(cat(2,resp_measure,' var explained'));
set(gca,'YLim',max(svd_s_dsq)*[1/logrange 1]);
%
subplot(2,3,2);
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
subplot(2,3,3);
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
axes('Position',[0.01,0.04,0.01,0.01]);
text(0,0,get(gcf,'Name'),'Interpreter','none');
axis off
%
if_done=0;
while (if_done==0)
    %
    %need to compute f-ratios for each, display frac of variance and f-ratios,
    %and choose pc;s based on pc#, frac var, f-rat, f-rat p-val
    %or set if_done=1
        npc_sel=getinp('number of pcs to select','d',[0 n_pc_max]);
        pc_sel=[1:npc_sel];
        pc_sel_string=sprintf('pc selection: up to %3.0f',npc_sel);
        if (npc_sel==0)
            if_done=1;
        end
    if (if_done==0)
        %
        proj_pc=svd_v(:,pc_sel)*svd_v(:,pc_sel)'; %projection on selected pcs
        vfilt_indiv_repts=v_indiv_repts*proj_pc;
        vfilt=reshape(reshape(v,[n_pixels*resp_maxlength,s.n_repts_kept*s.n_stims_kept])*proj_pc,size(v));
        %
        %make a heatmap of each slice, for each stimulus, averaged across time
        %
        if if_spatial_pattern
            v_perstim=reshape(mean(mean(reshape(vfilt_indiv_repts,[n_pixels,resp_minlength,s.n_repts_kept,s.n_stims_kept]),2),3),[s.n_pixels_kept,s.n_stims_kept]);
            v_perstim_range=[min(v_perstim(:)) max(v_perstim(:))];
            for plane_ptr=1:s.n_planes_with_data_kept
                plane=s.plane_list_kept(plane_ptr);
                pxls_inplane=find(s.xyz_kept(:,3)==plane);
                figure;
                set(gcf,'Position',[100 100 1200 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',cat(2,sprintf(' plane: %2.0f, ',plane),tstring));
                for stim_ptr_seq=1:s.n_stims_kept
                    stim_ptr=display_ptr_order(stim_ptr_seq);
                    stim_no=opts_read.stim_list(stim_ptr);
                    stim_name=stims.names_short{stim_no};
                    heatmap=NaN(x_len,y_len);
                    for ipxl=pxls_inplane'
                        heatmap(xyz_rel(ipxl,1),xyz_rel(ipxl,2))=v_perstim(ipxl,stim_ptr);
                    end
                    subplot(nrows,ncols,stim_ptr_seq);
                    imagesc(heatmap',v_perstim_range); %matlab plots the transpose
                    axis equal;
                    set(gca,'XTick',[1 x_len]);
                    set(gca,'XTickLabels',xyz_range(:,1));
                    set(gca,'YTick',[1 y_len]);
                    set(gca,'YTickLabels',xyz_range(:,2));
                    axis tight;
                    colorbar;
                    title(sprintf('s %2.0f: %s',stim_no,stim_name));
                end
                axes('Position',[0.01,0.01,0.01,0.01]);
                text(0,0,get(gcf,'Name'),'Interpreter','none');
                axis off
                axes('Position',[0.01,0.04,0.01,0.01]);
                text(0,0,pc_sel_string,'Interpreter','none');
                axis off
            end
        end
        %
        if if_temporal_pattern~=0
            v_pertime=reshape(mean(vfilt,1,'omitnan'),[resp_maxlength,s.n_repts_kept,s.n_stims_kept]); %average across space and repeat
            v_pertime_range=[min(v_pertime(:)) max(v_pertime(:))];
            %
            figure;
            set(gcf,'Position',[100 100 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'timecourse: ',tstring));
            for stim_ptr_seq=1:s.n_stims_kept
                stim_ptr=display_ptr_order(stim_ptr_seq);
                stim_no=opts_read.stim_list(stim_ptr);
                stim_name=stims.names_short{stim_no};
                subplot(nrows,ncols,stim_ptr_seq);
                plot(v_pertime(:,:,stim_ptr));
                hold on;
                plot(mean(v_pertime(:,:,stim_ptr),2,'omitnan'),'k','LineWidth',1);
                xlabel('response frame');
                set(gca,'XLim',[1 resp_maxlength]);
                set(gca,'YLim',v_pertime_range);
                title(sprintf('s %2.0f: %s',stim_no,stim_name));
            end
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,get(gcf,'Name'),'Interpreter','none');
            axis off
            axes('Position',[0.01,0.04,0.01,0.01]);
            text(0,0,pc_sel_string,'Interpreter','none');
            axis off
            if if_temporal_pattern==-1 %plot extended response with same vertical scale
                %
                %filter by pc
                v_alltime_recon=reshape(reshape(v_alltime,[s.n_timepoints_read*s.n_pixels_kept,s.n_repts_kept*s.n_stims_kept])*proj_pc,size(v_alltime));
                %
                figure;
                set(gcf,'Position',[100 100 1200 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',cat(2,'extended timecourse: ',tstring));
                for stim_ptr_seq=1:s.n_stims_kept
                    stim_ptr=display_ptr_order(stim_ptr_seq);
                    stim_no=opts_read.stim_list(stim_ptr);
                    stim_name=stims.names_short{stim_no};
                    %
                    v_stim=reshape(mean(v_alltime_recon(:,:,:,stim_ptr),1,'omitnan'),[s.n_timepoints_read s.n_repts_kept]);
                    %
                    subplot(nrows,ncols,stim_ptr_seq);
                    plot(v_stim);
                    hold on;
                    plot(mean(v_stim,2,'omitnan'),'k','LineWidth',1);
                    %plot stimulus onsets and offsets
                    colors=get(gca,'ColorOrder');
                    for ir=1:s.n_repts_kept
                        frame_on=s.onset_indexes_kept(ir,stim_ptr)/s.n_planes;
                        frame_off=s.offset_indexes_kept(ir,stim_ptr)/s.n_planes;
                        if if_timemark1
                            %method 1: cursors at each position
                            ypos=v_pertime_range(1)+(ir-1)*diff(v_pertime_range)/s.n_repts_kept+[0 diff(v_pertime_range)/s.n_repts_kept];
                            hp=plot(repmat(frame_on,1,2),ypos);
                            set(hp,'Color',colors(ir,:));
                            hp=plot(repmat(frame_off,1,2),ypos);
                            set(hp,'Color',colors(ir,:));
                        end
                        if if_timemark2
                            %method 2: bar for stim on
                            ypos=v_pertime_range(1)+(ir/4)*diff(v_pertime_range)/s.n_repts_kept;
                            hp=plot([frame_on frame_off],repmat(ypos,1,2),'k');
                            set(hp,'Color',colors(ir,:));
                        end
                    end
                    xlabel('frame');
                    set(gca,'XLim',[1 s.n_timepoints_read]);
                    set(gca,'YLim',v_pertime_range);
                    title(sprintf('s %2.0f: %s',stim_no,stim_name));
                end
                axes('Position',[0.01,0.01,0.01,0.01]);
                text(0,0,get(gcf,'Name'),'Interpreter','none');
                axis off
                axes('Position',[0.01,0.04,0.01,0.01]);
                text(0,0,pc_sel_string,'Interpreter','none');
                axis off
            end
        end
        %
        if if_dist
            %
            % average across repeats
            %
            v_across_repts=reshape(mean(vfilt,3,'omitnan'),[s.n_pixels_kept,resp_maxlength,s.n_stims_kept]);
            v_across_repts=reshape(v_across_repts(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_stims_kept]);
            %
            figure;
            set(gcf,'Position',[100 100 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'distances: ',tstring));
            for idist=1:n_dists
                switch dist_names{idist}
                    case 'cosine'
                        dot_prods=v_across_repts'*v_across_repts;
                        mags=sqrt(diag(dot_prods));
                        heatmap=dot_prods./(mags*mags');
                    case 'Pearson'
                        heatmap=corr(v_across_repts);
                end
                subplot(1,n_dists,idist);
                imagesc(heatmap(display_ptr_order,display_ptr_order)-diag(diag(heatmap))); %remove the diagonal so as not to inflate the scale
                set(gca,'XTick',[1:s.n_stims_kept]);
                set(gca,'XTickLabels',stims.names_short(s.opts_read.stim_list(display_ptr_order)));
                set(gca,'YTick',[1:s.n_stims_kept]);
                set(gca,'YTickLabels',stims.names_short(s.opts_read.stim_list(display_ptr_order)));
                axis square;
                title(dist_names{idist})
                colorbar;
            end
            clear v_across_repts
            %
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,cat(2,tstring,' ',rept_string),'Interpreter','none');
            axis off
            %
            % keep individual repeats
            %
            v_indiv_repts_dist=reshape(vfilt,[s.n_pixels_kept,resp_maxlength,s.n_repts_kept*s.n_stims_kept]);
            v_indiv_repts_dist=reshape(v_indiv_repts_dist(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_repts_kept*s.n_stims_kept]);
            %
            var_rats=hlid_varrats(reshape(v_indiv_repts_dist,[s.n_pixels_kept*resp_minlength,s.n_repts_kept,s.n_stims_kept]));
            %
            dpo_expanded=s.n_repts_kept*repmat(display_ptr_order-1,s.n_repts_kept,1); %expand the display pointer order to take into account blocks of repeats
            dpo_expanded=dpo_expanded+repmat([1:s.n_repts_kept]',1,s.n_stims_kept);
            dpo_expanded=dpo_expanded(:)';
            tick_locs=[1:s.n_stims_kept]*s.n_repts_kept-(s.n_repts_kept-1)/2; %tick locations for a block of repeats
            %
            figure;
            set(gcf,'Position',[100 100 1200 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'distances, each repeat: ',tstring));
            for idist=1:n_dists
                switch dist_names{idist}
                    case 'cosine'
                        dot_prods=v_indiv_repts_dist'*v_indiv_repts_dist;
                        mags=sqrt(diag(dot_prods));
                        heatmap=dot_prods./(mags*mags');
                    case 'Pearson'
                        heatmap=corr(v_indiv_repts_dist);
                end
                subplot(1,n_dists,idist);
                imagesc(heatmap(dpo_expanded,dpo_expanded)-diag(diag(heatmap))); %remove the diagonal so as not to inflate the scale
                set(gca,'XTick',tick_locs);
                set(gca,'XTickLabels',stims.names_short(s.opts_read.stim_list(display_ptr_order)));
                set(gca,'YTick',tick_locs);
                set(gca,'YTickLabels',stims.names_short(s.opts_read.stim_list(display_ptr_order)));
                axis square;
                title(dist_names{idist})
                colorbar;
            end
            %
            axes('Position',[0.01,0.07,0.01,0.01]);
            text(0,0,sprintf('F ratio: %8.4f',var_rats.frat));
            axis off
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,get(gcf,'Name'),'Interpreter','none');
            axis off
            axes('Position',[0.01,0.04,0.01,0.01]);
            text(0,0,pc_sel_string,'Interpreter','none');
            axis off
        end
    end
end %if_done


