%hlid_vi_pcafilt_coords_auto: create coordinates from KC volumetric imaging, data from George Barnum, Hong Lab,
%
% multiple methods of creating coordinates, from hlid_rastim_mds_coords_make.m
% derived from hlid_vi_pca_filt_auto,and can run over several datasets
% can customize pcrits
%  
%   See also:  HLID_VI_READ, HLID_VI_PCAFILT, HLID_VI_SPATIALFILTER, HLID_VI_STIMNAMES, HLID_VI_EXPLORE, HLID_VI_PCASELECT,
% HLID_METHS_DEFINE, HLID_RASTIM_MDS_COORDS_MAKE.
%
if_debug=getinp('1 for debug mode,','d',[0 1],0);
if ~exist('data_path') data_path='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\'; end
if ~exist('coord_file_infix') coord_file_infix='vi';end
if ~exist('meths')
    meths=hlid_meths_define;
end
if ~exist('n_stims') n_stims=24; end
if ~exist('n_repts') n_repts=5; end
if ~exist('max_timepoints') max_timepoints=0; end %until end of file
if ~exist('d_range_avg') d_range_avg=[-1 1]; end %colormap range for average responses
%
if_dist=1; % getinp('1 to plot distances','d',[0 1],1);
if_reorder=1; % getinp('1 to reorder stimuli','d',[0 1],1);
rept_list=[1:n_repts]; % getinp('repeat list','d',[1 n_repts],1:n_repts);
stim_list=[1:n_stims]; % getinp('stimulus list','d',[1 n_stims],1:n_stims);
%
if ~exist('opts_read') opts_read=struct;end
opts_read.max_timepoints=0; %read all time points
opts_read.if_remnan=1; %remove NaN's looking at all of the data
opts_read.data_path=data_path;
opts_read.stim_list=stim_list;
opts_read.rept_list=rept_list;
opts_read.if_keep_all_raw=0;
opts_read.if_log=0;
%
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
    '20250808_b_30s_output_fixed.hdf5',...
    '20250102_a_med_30s_output_mc_2.hdf5'};
end
data_files_med=find(contains(data_files,'_med'));
if_ok=0;
%starting values
data_files_selected=data_files_med;
if ~exist('sf_list') sf_list=[0:2]; end
resp_measures={'deltaF/F','z'};
if ~exist('resp_measure_list') resp_measure_list=1; end
if ~exist('meth_list') meth_list=[1 2 3 4]; end
%
if ~exist('pcrits') pcrits=[1 0.05]; end %critical values of p for pca filtering
 p_prefix='pr'; %for raw probability (could also be total prob, fdr corrected, etc)
ncols_input=1; %columns for pca plots
if ~exist('logrange') logrange=10^3; end %range of pc powers to plot
%
while (if_ok==0)
    for k=1:length(data_files)
        disp(sprintf(' %2.0f->%s',k,data_files{k}));
    end
    data_files_selected=getinp('choices','d',[1 length(data_files)],data_files_selected);
    n_files=length(data_files_selected);
    coord_file_infix=getinp('output file name infix','s',[],coord_file_infix);
    %
    sf_list=getinp('spatial filtering list full-widths','d',[0 6],sf_list);
    n_sfs=length(sf_list);
    %
    for k=1:length(resp_measures)
        disp(sprintf('%1.0f->response measure %s',k,resp_measures{k}));
    end
    resp_measure_list=getinp('choices','d',[1 length(resp_measures)],resp_measure_list);
    %
    pcrits=getinp('list of critical p-values for pca filtering','f',[0 1],pcrits);
    n_pcrits=length(pcrits);
    %
    for k=1:length(meths)
        disp(sprintf('%1.0f->%s',k,meths{k}.name_full));
    end
    meth_list=getinp('list of methods for dimension reduction','d',[1 length(meths)],meth_list);
    %
    coord_file_base=cell(1,n_files);
    for file_ptr=1:n_files
        ifile=data_files_selected(file_ptr);
        data_file=data_files{ifile};
        coord_file_id=cat(2,data_file(1:4),'-',data_file(5:6),'-',data_file(7:8),'-',data_file(10));
        coord_file_base{file_ptr}=cat(2,'hlid_',coord_file_infix,'_coords_',coord_file_id);
        disp(sprintf(' data file  %-40s will be used for coordinate files %s*.mat',data_file,coord_file_base{file_ptr}));
    end
    if_ok=getinp('1 if ok','d',[0 1]);
end
%
stims=hlid_vi_stimnames;
%
results=cell(n_files,n_sfs,length(resp_measure_list),n_pcrits,length(meth_list));
for file_ptr=1:n_files
    ifile=data_files_selected(file_ptr);
    data_file=data_files{ifile};
    %
    for sf_ptr=1:n_sfs
        sfilt=sf_list(sf_ptr);
        opts_read.if_spatialfilter=double(sfilt>0);
        opts_read.sfilt_hw=0.5*sfilt;
        opts_read.data_file=data_file;
        %
        disp('***********');
        disp(sprintf('processing file %3.0f of %3.0f: %-30s, reading with sf%1.0f',file_ptr,n_files,data_file,sfilt));
        %
        [s,opts_read_used]=hlid_vi_read(opts_read);
        n_pixels=s.n_pixels_kept;
        %
        read_data_file_short=strrep(s.opts_read.data_file,'.hdf5','');
        rept_string=cat(2,'repts: ',sprintf(' %2.0f',opts_read.rept_list));
        sf_tp_string=cat(2,sprintf('sf fw=%2.0f',sfilt),sprintf('; timepoints: [0 %2.0f]',s.n_timepoints_read));
        %
        tstring_input=cat(2,read_data_file_short,':',sf_tp_string,', ',rept_string);
        disp(sprintf('read %s',tstring_input));
        disp(sprintf('planes kept: %3.0f, total pixels kept: %7.0f, pixels per plane kept: %s',...
            s.n_pixels_kept,s.n_pixels_kept,sprintf(' %5.0f',s.pixels_per_plane_kept)))
        %
        if opts_read.if_log
            disp(s)
        end
        %
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
        for rm_ptr=1:length(resp_measure_list)
            rm=resp_measure_list(rm_ptr);
            resp_measure=resp_measures{rm};
            switch resp_measure
                case 'deltaF/F'
                    rm_string='dff';
                case 'z'
                    rm_string='z';
            end
            %
            %compute individual trial responses before pc filtering
            %
            %compute mean response measure (delta-F/F or z), averaged over repeats
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
            disp(sprintf('individual responses computed prior to pc filtering for response measure %s',resp_measure));
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
            tstring_input_resp=cat(2,tstring_input,' ',resp_measure);
            figh_pca=figure;
            set(gcf,'NumberTitle','off');
            set(gcf,'Position',[50 50 1400 800]);
            set(gcf,'Name',tstring_input_resp);
            %
            subplot(3,ncols_input,1)
            semilogy(svd_s_dsq,'.-');
            xlabel('eigenvalue');
            set(gca,'XLim',[0 n_pc_max]);
            ylabel(cat(2,resp_measure,' var explained'));
            set(gca,'YLim',max(svd_s_dsq)*[1/logrange 1]);
            %
            subplot(3,ncols_input,1+ncols_input);
            semilogy(frats);
            hold on;
            for icrit=1:n_pcrits
                finv_crit=finv(1-pcrits(icrit),varrats.fdof(1),varrats.fdof(2));
                semilogy([0 n_pc_max],repmat(finv_crit,1,2),'k:');
            end
            xlabel('eigenvalue');
            set(gca,'XLim',[0 n_pc_max]);
            ylabel('F ratio');
            set(gca,'YLim',[.1 100]);
            %
            subplot(3,ncols_input,1+2*ncols_input);
            semilogy(frats_pvals);
            hold on;
            for icrit=1:n_pcrits
                semilogy([0 n_pc_max],repmat(pcrits(icrit),1,2),'k:');
            end
            xlabel('eigenvalue');
            set(gca,'XLim',[0 n_pc_max]);
            ylabel('p(F ratio)');
            set(gca,'YLim',[10^-5 1]);
            %
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,tstring_input_resp,'Interpreter','none');
            axis off
            %
            %now analyze for each level of pca selection
            %
            figh_pcrits=figure;
            set(gcf,'NumberTitle','off');
            set(gcf,'Position',[50 50 1400 800]);
            set(gcf,'Name',cat(2,'pc filtering ',tstring_input_resp));
            [nrows_pcrits,ncols_pcrits]=nicesubp(n_pcrits,0.7);
            %
            stim_names=stims.names_short(s.opts_read.stim_list(display_ptr_order));
            for pcrit_ptr=1:length(pcrits)
                pcrit=pcrits(pcrit_ptr);
                %convert p=1 to all, p=0.05 to 050, p=0.001 to 001, p=.123 to 123, p=0.0001 to 00010
                if pcrit==1
                    pnum='all';
                elseif pcrit<0.001
                    pdec=5;
                    pnum=zpad(round(pcrit*10^pdec),pdec);
                else
                    pdec=3;
                    pnum=zpad(round(pcrit*10^pdec),pdec);
                end
                pcrit_string=cat(2,p_prefix,'-',pnum);
                %logic of hlid_vi_pcaselect, based on max p-value
                opts_pcasel=struct();
                opts_pcasel.n_pc_max=n_pc_max;
                opts_pcasel.eiv_squared=svd_s_dsq;
                opts_pcasel.frats=frats;
                opts_pcasel.fdof=varrats.fdof;
                pval_max=pcrit;
                frat_min=finv(1-pval_max,opts_pcasel.fdof(1),opts_pcasel.fdof(2));
                pc_sel=find(opts_pcasel.frats>=frat_min);
                frac_var_used=sum(opts_pcasel.eiv_squared(pc_sel))/sum(opts_pcasel.eiv_squared(:));
                pc_sel_string=sprintf('p %5.3f F %5.3f npc %3.0f pfrac %5.3f',pval_max,frat_min,length(pc_sel),frac_var_used);
                disp(sprintf('pca set %2.0f: %s',pcrit_ptr,pc_sel_string));
                %
                %reconstruct from selected pc;s
                proj_pc=svd_v(:,pc_sel)*svd_v(:,pc_sel)'; %projection on selected pcs
                vfilt_indiv_repts=v_indiv_repts*proj_pc;
                vfilt=reshape(reshape(v,[n_pixels*resp_maxlength,s.n_repts_kept*s.n_stims_kept])*proj_pc,size(v));
                % average across repeats
                v_across_repts=reshape(mean(vfilt,3,'omitnan'),[s.n_pixels_kept,resp_maxlength,s.n_stims_kept]);
                v_across_repts=reshape(v_across_repts(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_stims_kept]);
                % indiv repeats to compute F
                v_indiv_repts_dist=reshape(vfilt,[s.n_pixels_kept,resp_maxlength,s.n_repts_kept*s.n_stims_kept]);
                v_indiv_repts_dist=reshape(v_indiv_repts_dist(:,[1:resp_minlength],:),[s.n_pixels_kept*resp_minlength,s.n_repts_kept*s.n_stims_kept]);
                varrats_indiv=hlid_varrats(reshape(v_indiv_repts_dist,[s.n_pixels_kept*resp_minlength,s.n_repts_kept,s.n_stims_kept]));
                %
                dot_prods=v_across_repts'*v_across_repts;
                mags=sqrt(diag(dot_prods));
                heatmap=dot_prods./(mags*mags');
                figure(figh_pcrits);
                subplot(nrows_pcrits,ncols_pcrits,pcrit_ptr);
                heatmap_ordered=heatmap(display_ptr_order,display_ptr_order);
                imagesc(heatmap_ordered-diag(diag(heatmap_ordered)),d_range_avg); %remove the diagonal so as not to inflate the scale
                colorbar;
                set(gca,'XTick',[1:s.n_stims_kept]);
                set(gca,'XTickLabels',stim_names);
                set(gca,'YTick',[1:s.n_stims_kept]);
                set(gca,'YTickLabels',stim_names);
                xlabel(sprintf('overall f-ratio: %6.3f',varrats_indiv.frat));
                axis square;
                title(pc_sel_string,'FontSize',8,'FontWeight','normal');
                %
                for meth_ptr=1:length(meth_list)
                    meth=meth_list(meth_ptr);
                    meth_string=meths{meth}.name_file;
                    %
                    coord_file=cat(2,coord_file_base{file_ptr},'_sf',sprintf('%1.0f',sfilt),'_',rm_string,'_',pcrit_string,'_',meth_string);
                    if if_debug>0
                        disp(coord_file)                       
                    end
                    r=struct;
                    r.data_file=data_file;
                    r.coord_file=coord_file;
                    r.stims=stims;
                    %
                    r.sfilt=sfilt;
                    r.resp_measure=resp_measure;
                    r.eiv_squared=svd_s_dsq;
                    %
                    r.pcrit=pcrit;
                    r.pc_sel=pc_sel;
                    r.frat_min=frat_min;%min f-ratio for a component to be kept
                    r.frac_var_used=frac_var_used;
                    r.frat_overall=varrats_indiv.frat;
                    %
                    r.meth_string=meth_string;
                    %
                    results{file_ptr,sf_ptr,rm_ptr,pcrit_ptr,meth_ptr}=r;
                end %dim reduction method
            end %pcrit_ptr
            figure(figh_pcrits);
            axes('Position',[0.01,0.01,0.01,0.01]);
            text(0,0,cat(2,'dot products ',tstring_input_resp),'Interpreter','none');
            axis off
        end %rm_ptr (response measure)
        clear s
    end % sf_ptr
end %file_ptr
