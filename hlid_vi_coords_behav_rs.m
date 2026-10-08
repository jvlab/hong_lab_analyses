% hlid_vi_coords_behav_rs: read volumetric imaging coordinates set files, compare with behavioral data
%
% behavioral data (valence) indicated by color
% reliability of behavior across flies indicated by size (large: less variable)
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, RS_EXTRACT_COORDSETS, HLID_VI_COORDS_KNIT_RS, HLID_VI_COORDS_KNIT_RS_AUTO, REGRESS.
%
hlid_setup;
%
if ~exist('path_beh') path_beh='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC'; end
if ~exist('filename_beh') filename_beh='attractiveness_per_fly_overview.xlsx'; end
%
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
table_beh=readtable(cat(2,path_beh,filesep,filename_beh));
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
%
typenames_beh=table_beh.Properties.VariableDescriptions(2:end)'; %skip the index column
nstims_beh=length(typenames_beh);
rownames_beh=table_beh.fly;
rownums_beh=strmatch('202',rownames_beh); %fly names are 2024_*, etc
%
data_beh=table2array(table_beh(rownums_beh,2:end))'; %rows are stimuli, cols are preps
nstims_beh=size(data_beh,1);
npreps_beh=size(data_beh,2);
disp(sprintf('behavioral data read for %2.0f stimuli from %2.0f preps',nstims_beh,npreps_beh));
%
if ~exist('tol') tol=10^-5; end %for matching stat calcs with matlab
if ~exist('dmax') dmax=5; end
if ~exist('dims_plot') dims_plot=[2 3];end
if ~exist('line_width') line_width=2; end
if ~exist('marker_size') marker_size=12; end
%
dmax=getinp('maximum dimension to analyze','d',[2 10],dmax);
dims_plot=getinp('dimensions to plot (0 for none)','d',[0 dmax],dims_plot);
dims_plot=setdiff(dims_plot,[0 1]);
if ~isempty(dims_plot)
    if_plain=getinp('1 to also plot spaces without color-by-behavior','d',[0 1],0);
end
%
if ~exist('opts_read') opts_read=struct(); end
opts_read.input_type=1; %just data
opts_read.if_warn=0; %may have different sets of stimuli or stimuli in different orders
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.ui_filter='*coords*consensus*sf0*dff*pr-all*';
opts_read.if_log=0;
%
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
%
if ~exist('opts_disp') opts_disp=struct(); end
opts_disp.coord_group_method='keeplow';
opts_disp.if_legend=0;
opts_disp.callout_amount=0.5;
opts_disp.set_markersizes=marker_size;
opts_disp.set_colors='k';
%
aux=struct;
aux.nsets=0;
aux.opts_check=opts_check;
aux.opts_read=opts_read;
%
fullnames=[];
[data_read,aux_read_out]=rs_get_coordsets(fullnames,aux);
nsets=length(data_read.sets);
disp(sprintf(' %3.0f sets read.',nsets));
mean_beh=mean(data_beh,2);
range_mean_beh=[min(mean_beh),max(mean_beh)];
if ~exist('colors_beh') colors_beh=[1 0 0;0 1 1]; end
%
std_beh=std(data_beh,0,2);
range_std_beh=[min(std_beh),max(std_beh)];
if ~exist('sizes_beh') sizes_beh=[30;15]; end %render std dev by symbol size (smaller std is bigger symbol)
%
for iset=1:nsets
    typenames=data_read.sas{iset}.typenames;
    nstims=length(typenames);
    ptrs_beh=zeros(1,nstims);
    for k=1:nstims_beh
        ptr_beh=strmatch(typenames_beh{k},typenames,'exact');
        if length(ptr_beh)~=1
            ptr_beh=0;
        end
        ptrs_beh(k)=ptr_beh;
    end
    filename=strrep(strrep(data_read.sets{iset}.label_long,'/',filesep),'\',filesep);
    filename_end=max(find(cat(2,filesep,filename)==filesep));
    filename_short=filename(filename_end:end);
    filename_short=strrep(filename_short,'.mat','');
    disp(sprintf('analyzing %s',filename_short));   
    if any(ptrs_beh==0)
        disp(sprintf('no coords for %1.0f behaviors found in %s',sum(ptrs_beh==0),filename_short));
    end
    y=mean_beh(find(ptrs_beh>0)); %regress agains behaviors that have coordinates
    pvals=zeros(dmax,1);
    frats=zeros(dmax,1);
    rsquareds=zeros(dmax,1);
    rsquareds_drop=zeros(dmax,1);
    rmse=zeros(dmax,1);
    rmse_drop=zeros(dmax,1);
    b=cell(dmax,1);
    b_intvl=cell(dmax,1);
    b_drop=cell(dmax,1);
    sig_string=cell(dmax,1);
    for dim=1:dmax
        coords=data_read.ds{iset}{dim};
        x=coords(ptrs_beh(ptrs_beh>0),:); %regress against behaviors that have coords
        n=size(x,1);
        [b{dim},b_intvl{dim},r,r_intvl,stats]=regress(y,[ones(n,1),x]); %add a constant term
        disp(sprintf(' dim %1.0f: regressors (constant and each pc), and 0.95 confidence limits',dim))
        disp([b{dim},b_intvl{dim}]')
        for k=1:dim
            if sign(b_intvl{dim}(k+1,1))==sign(b_intvl{dim}(k+1,2))
                sig_string{dim}=cat(2,sig_string{dim},sprintf(' dim %2.0f ',k));
            else
                sig_string{dim}=cat(2,sig_string{dim},'        ');
            end
        end
        %stats: the R-square statistic, the F statistic, p value for the full model, and an estimate of the error variance.
%       disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, from stats',stats(3),stats(2),stats(1)));
        %recalculate stats, first principles
        y_pred=[ones(n,1),x]*b{dim};
        ss_model=sum((y_pred-mean(y)).^2);
        ss_error=sum((y-y_pred).^2);
        frats(dim)=(ss_model/dim)/(ss_error/(n-dim-1));
        pvals(dim)=1-fcdf(frats(dim),dim,n-dim-1);
        rsquareds(dim)=corr(y,y_pred).^2;
%       disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, recalc',p,frat,Rsquared));
        %
        if abs(pvals(dim)-stats(3))>tol
            disp(sprintf('mismatch of p: %7.3f (matlab) vs %7.3f (recalc)',pvals(dim),stats(3)));
        end
        if abs(frats(dim)-stats(2))>tol
            disp(sprintf('mismatch of F-ratio: %7.3f (matlab) vs %7.3f (recalc)',frats(dim),stats(2)));
        end
        if abs(rsquareds(dim)-stats(1))>tol
            disp(sprintf('mismatch of R-squared: %7.3f (matlab) vs %7.3f (recalc)',rsquareds(dim),stats(1)));
        end
        y_pred_drop=zeros(n,1);
        for k=1:n
            i_drop=setdiff([1:n],n);
            x_drop=x(i_drop,:);
            y_drop=y(i_drop,:); 
            b_drop{dim}=regress(y_drop,[ones(n-1,1),x_drop]);
            y_pred_drop(k)=[1 x(k,:)]*b_drop{dim};
        end
        ss_error_drop=sum((y-y_pred_drop).^2);
        rmse(dim)=sqrt(ss_error/n);
        rmse_drop(dim)=sqrt(ss_error_drop/n);
        rsquareds_drop(dim)=corr(y,y_pred_drop).^2;
    end %dim
    disp('   dim     p    f-ratio      R^2    rmse     R^2_drop rmse_drop   signif regressors');
    for dim=1:dmax
        disp(sprintf('%5.0f  %7.3f %7.4f    %7.4f %7.4f    %7.4f %7.4f         %s',dim,pvals(dim),frats(dim),rsquareds(dim),rmse(dim),rsquareds_drop(dim),rmse_drop(dim),sig_string{dim}))
    end
    std_beh_have=std(data_beh(ptrs_beh>0,:),0,2);
    rms_beh=sqrt(mean(std_beh_have.^2));
    disp(sprintf(' rms dev for behavior: %7.4f',rms_beh));
    %
    %plot
    %
    if ~isempty(dims_plot)
        for dim_plot=dims_plot
            opts_disp.set_select=1;
            opts_disp.dim_select=dim_plot;
            %
            if if_plain
                %simple plot, all points black
                data_read_oneset=rs_extract_coordsets(data_read,iset);
                aux_disp=rs_disp_coordsets(data_read_oneset,setfield(aux,'opts_disp',opts_disp));
                set(gcf,'Name',sprintf('dim %1.0f, %s',dim_plot,filename_short));
                %
                axes('Position',[0.01,0.04,0.01,0.01]); %for text
                text(0,0,filename_short,'Interpreter','none','FontSize',8);
                axis off;
            end
            %
            %plot with custom colors for each point: make each point a different set
            %callouts will have slightly different lengths, since they are normalized by rms within each set
            %
            ptrs_beh_have=ptrs_beh(ptrs_beh>0);
            %
            data_indiv=struct;
            data_indiv.sets=cell(1,n);
            data_indiv.sas=cell(1,n);
            data_indiv.ds=cell(1,n);
            for istim_beh=1:n
                istim=ptrs_beh_have(istim_beh); %which stimulus, in list of typenames
                data_indiv.sets{istim_beh}=data_read.sets{iset};
                data_indiv.sets{istim_beh}.nstims=1;
                data_indiv.sas{istim_beh}=data_read.sas{iset};
                data_indiv.sas{istim_beh}.nstims=1;
                data_indiv.sas{istim_beh}.typenames=data_read.sas{iset}.typenames(istim);
                data_indiv.sas{istim_beh}.btc_specoords=data_read.sas{iset}.btc_specoords(istim,:);
                data_indiv.ds{istim_beh}=data_read.ds{iset};
                istim=ptrs_beh_have(istim_beh);
                for dq=1:length(data_indiv.ds{istim_beh})
                    data_indiv.ds{istim_beh}{dq}=data_indiv.ds{istim_beh}{dq}(istim,:);
                end
                % istim_beh
                % istim
                % typenames(istim)
            end
            opts_disp_indiv=opts_disp;
            opts_disp_indiv.data_label_setsel_method='all';
            opts_disp_indiv.set_select=[1:n]; %each set is a different stimulus
            opts_disp_indiv.data_label_typenames_vary=1;
            opts_disp_indiv.callout_center=mean(data_read.ds{iset}{dim_plot},1);
            opts_disp_indiv.callout_colors='set_colors';
            opts_disp_indiv.set_colors=cell(1,n);
            %assign behavior value to each point based on data_indiv.sas{istim}.typenames          
            for k=1:n
                frac_beh=(y(k)-range_mean_beh(1))/diff(range_mean_beh);
                opts_disp_indiv.set_colors{k}=[(1-frac_beh) frac_beh]*colors_beh;
                % [y(index_beh) opts_disp_indiv.set_colors{k}]
                istim=ptrs_beh_have(k);
                frac_std=(std_beh(istim)-range_std_beh(1))/diff(range_std_beh);
                size_std=[(1-frac_std) frac_std]*sizes_beh;
                opts_disp_indiv.set_markersizes(k)=round(size_std);
                % [k std_beh(istim) frac_std size_std]
                % typenames{k}
            end
            %
            aux_disp_indiv=rs_disp_coordsets(data_indiv,setfield(aux,'opts_disp',opts_disp_indiv));
            set(gcf,'Name',sprintf('dim %1.0f, %s',dim_plot,filename_short));
            %
            %plot vectors corresponding to what the regressors project on
            %
            coord_groups=aux_disp_indiv.opts_disp.coord_groups;
            axis_handles=aux_disp_indiv.opts_disp.axis_handles;
            for isub=1:length(axis_handles)
                axes(axis_handles{isub});
                hold on;
                cg=coord_groups(isub,:);
                vec=b{dim_plot}(1+cg); %first entry of b is constant term
                coords=data_read.ds{iset}{dim_plot}(:,cg);
                coord_mean=mean(coords,1);
                coord_scale=sqrt(mean(coords(:).^2)); %scale of coordinates
                vec_scale=coord_scale*vec./sqrt(sum(vec.^2)); %normalized
                for im=1:2
                    im_sign=-3+2*im;
                    if dim_plot>=3
                        hp=plot3(coord_mean(1)+[0 im_sign*vec_scale(1)],coord_mean(2)+[0 im_sign*vec_scale(2)],coord_mean(3)+[0 im_sign*vec_scale(3)],'k');
                    else
                        hp=plot(coord_mean(1)+[0 im_sign*vec_scale(1)],coord_mean(2)+[0 im_sign*vec_scale(2)],'k');
                    end
                    set(hp,'LineWidth',line_width);
                    set(hp,'Color',colors_beh(im,:));
                end
            end
            %
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,filename_short,'Interpreter','none','FontSize',8);
            axis off;
        end %dims_plot
    end %dims_plot empty?
end %iset
