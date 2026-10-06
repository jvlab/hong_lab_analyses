% hlid_vi_coords_behav_rs: read volumetric imaging coordinates set files, compare with behavioral data
%
%to do: graphics -- show regression direction in rep space
% 
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_COORDS_KNIT_RS, HLID_VI_COORDS_KNIT_RS_AUTO, REGRESS.
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
if ~exist('dims_plot') dims_plot=[3];end
dmax=getinp('maximum dimension to analyze','d',[2 10],dmax);
dims_plot=getinp('dimensions to plot (0 for none)','d',[0 7],dims_plot);
dims_plot=setdiff(dims_plot,[0 1]);
%
if ~exist('opts_read') opts_read=struct(); end
opts_read.input_type=1; %just data
opts_read.if_warn=0; %may have differnt sets of stimuli
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
opts_disp.set_markersizes=12;
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
range_beh=[min(mean_beh),max(mean_beh)];
if ~exist('colors_beh') colors_beh=[1 0 0;0 1 1]; end
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
    for dim=1:dmax
        coords=data_read.ds{iset}{dim};
        x=coords(ptrs_beh(ptrs_beh>0),:); %regress against behaviors that have coords
        n=size(x,1);
        [b,b_intvl,r,r_intvl,stats]=regress(y,[ones(n,1),x]); %add a constant term
        disp(sprintf(' dim %1.0f: regressors (constant and each pc), and 0.95 confidence limits',dim))
        disp([b,b_intvl]')
        %stats: the R-square statistic, the F statistic, p value for the full model, and an estimate of the error variance.
%       disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, from stats',stats(3),stats(2),stats(1)));
        %recalculate stats, first principles
        y_pred=[ones(n,1),x]*b;
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
            b_drop=regress(y_drop,[ones(n-1,1),x_drop]);
            y_pred_drop(k)=[1 x(k,:)]*b_drop;
        end
        ss_error_drop=sum((y-y_pred_drop).^2);
        rmse(dim)=sqrt(ss_error/n);
        rmse_drop(dim)=sqrt(ss_error_drop/n);
        rsquareds_drop(dim)=corr(y,y_pred_drop).^2;
    end
    disp('   dim     p    f-ratio      R^2    rmse     R^2_drop rmse_drop');
    for dim=1:dmax
        disp(sprintf('%5.0f  %7.3f %7.4f    %7.4f %7.4f    %7.4f %7.4f',dim,pvals(dim),frats(dim),rsquareds(dim),rmse(dim),rsquareds_drop(dim),rmse_drop(dim)))
    end
    std_beh=std(data_beh(ptrs_beh>0,:),0,2);
    rms_beh=sqrt(mean(std_beh.^2));
    disp(sprintf(' rms dev for behavior: %7.4f',rms_beh));
    %
    %plot
    %
    if ~isempty(dims_plot)
        opts_disp.set_select=iset;
        for dim_plot=dims_plot
            opts_disp.dim_select=dim_plot;
            aux_disp=rs_disp_coordsets(data_read,setfield(aux,'opts_disp',opts_disp));
            %
            %make each point a different set
            data_indiv=struct;
            data_indiv.sets=cell(1,n);
            data_indiv.sas=cell(1,n);
            data_indiv.ds=cell(1,n);
            for istim=1:n
                data_indiv.sets{istim}=data_read.sets{iset};
                data_indiv.sets{istim}.nstims=1;
                data_indiv.sas{istim}=data_read.sas{iset};
                data_indiv.sas{istim}.nstims=1;
                data_indiv.sas{istim}.typenames=data_read.sas{iset}.typenames(istim);
                data_indiv.sas{istim}.btc_specoords=data_read.sas{iset}.btc_specoords(istim,:);
                data_indiv.ds{istim}=data_read.ds{iset};
                for dq=1:length(data_indiv.ds{istim})
                    data_indiv.ds{istim}{dq}=data_indiv.ds{istim}{dq}(istim,:);
                end
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
                index_beh=strmatch(data_indiv.sas{k}.typenames{1},typenames_beh,'exact');
                frac_beh=y(index_beh)-range_beh(1)/diff(range_beh);
                opts_disp_indiv.set_colors{k}=[(1-frac_beh) frac_beh]*colors_beh;
                % [y(index_beh) opts_disp_indiv.set_colors{k}]
            end
            %
            aux_disp_indiv=rs_disp_coordsets(data_indiv,setfield(aux,'opts_disp',opts_disp_indiv));
            %
            coord_groups=aux_disp.opts_disp.coord_groups; %will need this to show regression vectors
            %
            axes('Position',[0.01,0.04,0.01,0.01]); %for text
            text(0,0,filename_short,'Interpreter','none','FontSize',8);
            axis off;
        end
    end
end
