%hlid_majaxes_compare: compare major axes of a transformatoin and axes determined from 
% connectivity matrix
% 
% uses a file saved by hlid_majaxes, containing axes of a transformation
% for plots in stimulus space, kc order is used
%
%   See also: HLID_MAJAXES, HLID_CSV2CONNECTIVTY_DEMO, PSG_MAJAXES.
%
hlid_setup;
plot_order=display_orders.kcclust; %stimulus plot order
model_types_def=psg_geomodels_define();
%
if ~exist('fn_majaxes') fn_majaxes='hlid_majaxes_roi_orn9setsMerge_kcConsensusScale_d6.mat'; end
if ~exist('fn_connectivity') fn_connectivity='hlid_csv2connectivity_19Sep24.mat'; end
if ~exist('neivs_connect') neivs_connect=4; end %number of eigenvecs of connectivity to show
if ~exist('nshuffs') nshuffs=1000; end %number of shuffles to compute
if ~exist('maxcorr_plot') maxcorr_plot=0.75; end %max correlation to plot
if ~exist('quantiles_plot') quantiles_plot=[0.5 0.05 0.01]; end
if ~exist('if_frozen') if_frozen=1; end %frozen random number gens
%
fn_connectivity=getinp('file name with connectivity information','s',[],fn_connectivity);
load(fn_connectivity);
%
connectivity_fields=fieldnames(connectivity);
c_filenames=cell(0);
c_varnames=cell(0);
c_fill_methods=cell(0);
for ic=1:length(connectivity_fields)
    icf=connectivity_fields{ic};
    c_filenames{end+1}=connectivity.(icf).filename;
    c_varnames{end+1}=connectivity.(icf).v_name;
    c_fill_methods{end+1}=connectivity.(icf).fill_method;
end
%
c_filenames=unique(c_filenames);
c_varnames=unique(c_varnames);
c_fill_methods=unique(c_fill_methods);
%
fn_majaxes=getinp('file name with axes information from transformations (e.g., hlid_majaxes_roi*.mat)','s',[],fn_majaxes);
load(fn_majaxes)
%code borrowed from psg_majaxes
nds_ref=size(results_axes,1);
nds_adj=size(results_axes,2);
d_ref_list=[];
d_adj_list=[];
%
[coords_adj,sa_adj]=psg_read_coorddata(fn_adj,[],opts_read); %also read original coordinates
%
for id_ref=1:nds_ref
    for id_adj=1:nds_adj
        res=results_axes{id_ref,id_adj};
        if ~isempty(res)
            d_ref_list=[d_ref_list,res.ref_dim];
            d_adj_list=[d_adj_list,res.adj_dim];
            model_types=res.model_types;
            ref_file=res.ref_file;
            adj_file=res.adj_file;
            results{id_ref,id_adj}.ref_dim=res.ref_dim;
            results{id_ref,id_adj}.adj_dim=res.adj_dim;
            results{id_ref,id_adj}.ref_file=ref_file;
            results{id_ref,id_adj}.adj_file=adj_file;
        end
    end %id_adj
end %id_ref
d_ref_list=unique(d_ref_list);
d_adj_list=unique(d_adj_list);
disp(sprintf('reference dataset: %s',ref_file));
disp('dimensions available:');
disp(d_ref_list);
disp(sprintf('adjusted dataset: %s',adj_file));
disp('dimensions available:');
disp(d_adj_list);
disp('model types used for fitting geometric transformation from adjusted to reference');
disp(model_types);
%calculations
adj_ref_labels={'adj','ref'};
nar=length(adj_ref_labels);
%
if_done=0;
while (if_done==0)
    % specify a transformation: dimension in adjusted space, reference space, model type,and component
    d_adj=getinp('adjusted dimension','d',[min(d_adj_list),max(d_adj_list)]);
    d_adj_ptr=find(d_adj_list==d_adj);
    d_ref=getinp('reference dimension','d',[min(d_ref_list),max(d_ref_list)]);
    d_ref_ptr=find(d_ref_list==d_ref);
    if ~isempty(d_adj_ptr) & ~isempty(d_ref_ptr)
        for im_ptr=1:length(model_types)
            disp(sprintf('%1.0f->%s',im_ptr,model_types{im_ptr}));
        end
        im_ptr=getinp('model type','d',[1 length(model_types)],1);
        trans_data=results_axes{d_ref_ptr,d_adj_ptr}.adj;
        projections_roi_all=trans_data.projections_roi{im_ptr};
        roi_names=trans_data.roi_names;
        ipw=getinp('piecewise component','d',[1 size(projections_roi_all,3)],1);
        projections_roi=projections_roi_all(:,:,ipw); %these are the directions identifed in the transformation
        projections_label=sprintf('%s: adj %3.0f ref %2.0f ipw %2.0f',model_types{im_ptr},d_adj,d_ref,ipw);
        %
        magnifs=trans_data.magnifs{im_ptr}(:,ipw);
        projections=trans_data.projections{im_ptr}(:,:,ipw);
        ytl_mags=cell(d_adj,1);
        for k=1:d_adj
            ytl_mags{k}=sprintf('dir %1.0f: %5.2f',k,magnifs(k));
        end
        %specify connectivity data
        for k=1:length(c_filenames)
            disp(sprintf('connectivity file %1.0f: %s',k,c_filenames{k}));
        end
        c_filename=c_filenames{getinp('choice','d',[1 length(c_filenames)])};
        for k=1:length(c_varnames)
            disp(sprintf('variable %1.0f: %s',k,c_varnames{k}));
        end
        c_varname=c_varnames{getinp('choice','d',[1 length(c_varnames)])};
        for k=1:length(c_fill_methods)
            disp(sprintf('fill method %1.0f: %s',k,c_fill_methods{k}));
        end
        c_fill_method=c_fill_methods{getinp('choice','d',[1 length(c_fill_methods)])};
        c_desc=sprintf('file: %s var: %s, fill method: %s',c_filename,c_varname,c_fill_method);
        c_string=cat(2,c_filename,'__',c_varname,'__',c_fill_method);
        if isfield(connectivity,c_string)
            c_data=connectivity.(c_string);
        else
            disp('not found in connectivity file')
        end
        %
        %match the roi_names with connectivity names
        %and vc3 matches with sum of vc3l and vc3m
        %
        if_roimatch=1;
        conn_use=zeros(length(roi_names),length(c_data.rois));
        for k=1:length(roi_names)
            roi_match=strmatch(roi_names{k},c_data.rois,'exact');
            if length(roi_match)==0
                if strcmp(roi_names{k},'VC3')
                    roi_matches=strmatch(roi_names{k},c_data.rois); %allow VC3l and VC3m to match
                    conn_use(k,roi_matches)=1;
                else
                    disp(sprintf('no match for roi %s',roi_names{k}));
                    if_roimatch=0;
                end
            elseif length(roi_match)==1
                conn_use(k,roi_match)=1;
            else
                disp(sprintf('multiple matches for roi %s',roi_names{k}));
                if_roimatch=0;
            end
        end
        %
        %compute correlations with projections_roi and shuffles
        %
        if if_frozen
            rng('default');
        end
        if if_roimatch
           conn_eivecs=conn_use*c_data.eigenvectors(:,1:neivs_connect); %eigenvectors of connectivity matrix, mapped to roi
            corrs_data=abs(corr(projections_roi,conn_eivecs));
            corrs_shuff=zeros(d_adj,neivs_connect,nshuffs);
            for k=1:nshuffs
                corrs_shuff(:,:,k)=abs(corr(projections_roi(randperm(size(projections_roi,1)),:),conn_eivecs));
            end
        end
        %
        figure;
        figure_name=cat(2,'connectivity: ',c_string,' transform: ',projections_label);
        set(gcf,'Position',[50 50 1400 900]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',figure_name);
        %
        ncols=4;
        nrows=3;
        %plot original coordinates
        %subtracting mean, expanding to full scale, and reordering
        subplot(nrows,ncols,1);
        z=coords_adj{d_adj};
        z=z-repmat(mean(z,1),size(z,1),1);
        [zplot,xtick_labels]=psg_majaxes_reorder(z,plot_order,trans_data.typenames);
        imagesc(zplot',max(abs(z(:)))*[-1 1]);
        ytl_coords=cell(1,d_adj);
        for k=1:d_adj
            ytl_coords{k}=sprintf('coord %1.0f: rms %5.2f',k,sqrt(sum(z(:,k).^2))); %plot rms 
        end
        set(gca,'XTick',[1:length(trans_data.typenames)]);
        set(gca,'XTickLabel',xtick_labels);
        set(gca,'YTick',[1:d_adj]);
        set(gca,'YTickLabel',ytl_coords);
        title('adj coords')        
        %plot projections onto odor space
        %subtracting mean, expanding to full scale, and reordering 
        subplot(nrows,ncols,2);
        z=projections;
        z=z-repmat(mean(z,1),size(z,1),1);
        [zplot,xtick_labels]=psg_majaxes_reorder(z,plot_order,trans_data.typenames);
        imagesc(zplot',max(abs(z(:)))*[-1 1]);
        set(gca,'XTick',[1:length(trans_data.typenames)]);
        set(gca,'XTickLabel',xtick_labels);
        set(gca,'YTick',[1:d_adj]);
        set(gca,'YTickLabel',ytl_mags);
        title('odor space xform dirs')
        %roi space
        subplot(nrows,ncols/2,2)
        z=projections_roi;
        z=z-repmat(mean(z,1),size(z,1),1);
        imagesc(z',max(abs(z(:)))*[-1 1]);
        set(gca,'XTick',[1:length(roi_names)]);
        set(gca,'XTickLabel',roi_names);
        set(gca,'YTick',[1:d_adj]);
        set(gca,'YTickLabel',ytl_mags);
        title('roi space xform dirs')
        %
        % show connectivity matrix filled in
        %
        subplot(nrows,ncols,ncols+1);
        imagesc(c_data.vals_filled);
        axis square;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        title('connectivity filled')
        if if_roimatch %further calculations if all rois match to connection tags
            %
            % show connectivity selection and reordering if all rois match up
            %
            subplot(nrows,ncols,ncols+2)
            imagesc(conn_use*c_data.vals_filled*conn_use');
            axis square;
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            title('connectivity at rois')
            %
            %show eigenvectors of connectivity matrix
            %
            z=conn_eivecs;
            subplot(nrows,ncols/2,ncols/2+2);
            z=z-repmat(mean(z,1),size(z,1),1);
            imagesc(z',max(abs(z(:)))*[-1 1]);
            ytl_ceiv=cell(1,neivs_connect);
            for k=1:neivs_connect
                ytl_ceiv{k}=sprintf('conn %1.0f: %5.2f',k,c_data.eigenvalues(k));
            end
            set(gca,'XTick',[1:length(roi_names)]);
            set(gca,'XTickLabel',roi_names);
            set(gca,'YTick',[1:neivs_connect]);
            set(gca,'YTickLabel',ytl_ceiv);
            title('connectivity eigenvectors');
            %
            %show correlations of transform dirs and connectivity dirs
            %
            for k=1:d_adj
                subplot(nrows,max(ncols,d_adj),k+max(ncols,d_adj)*(nrows-1));
                legt=cell(1,1+length(quantiles_plot));
                plot([1:neivs_connect],corrs_data(k,:),'k');
                legt{1}='data';
                hold on;
                for iq=1:length(quantiles_plot)
                    plot([1:neivs_connect],quantile(corrs_shuff(k,:,:),1-quantiles_plot(iq),3),'k--');
                    legt{iq+1}=sprintf('q %5.3f',quantiles_plot(iq));
                end
                set(gca,'XTick',[1:neivs_connect]);
                hold on;
                xlabel('conn eivec');
                ylabel('correl');
                set(gca,'XLim',0.5+[0 neivs_connect]);
                set(gca,'YLim',[0 max(maxcorr_plot,max(corrs_data(:)))]);
                legend(legt,'FontSize',7,'Location','Best');
                title(sprintf('xform dir %1.0f',k));
            end
        end
        %
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,sprintf('connectivity data from %s: %s  (%s)',fn_connectivity,c_string,c_desc),'Interpreter','none');
        axis off;
        axes('Position',[0.01,0.04,0.01,0.01]); %for text
        text(0,0,sprintf('major transformation axes from %s: %s',fn_majaxes,projections_label),'Interpreter','none');
        axis off;
        axes('Position',[0.01,0.06,0.01,0.01]); %for text
        text(0,0,sprintf('nshuffs: %7.0f',nshuffs));
        axis off;
        %
        colormap jet;
    else
        disp('not found in majaxes file')
    end
    if_done=getinp('1 if done','d',[0 1]);
end
