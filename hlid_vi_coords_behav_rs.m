% hlid_vi_coords_behav_rs: read volumetric imaging coordinates set files, compare with behavioral data
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
if ~exist('opts_read') opts_read=struct(); end
opts_read.input_type=1; %just data
opts_read.if_warn=0; %may have differnt sets of stimuli
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.ui_filter='*coords*consensus*sf0*dff*pr-all*';
opts_read.if_log=0;
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
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
    for dim=1:5
        coords=data_read.ds{iset}{dim};
        x=coords(ptrs_beh(ptrs_beh>0),:); %regress against behaviors that have coords
        n=size(x,1);
        [b,b_intvl,r,r_intvl,stats]=regress(y,[ones(n,1),x]); %add a constant term
        disp(sprintf(' dim %1.0f: regrssors (constant and each pc), and 95% confidence limits',dim))
        disp([b,b_intvl]')
        %stats: the R-square statistic, the F statistic, p value for the full model, and an estimate of the error variance.
        disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, from stats',stats(3),stats(2),stats(1)));
        %recalculate stats, first principles
        y_pred=[ones(n,1),x]*b;
        ss_model=sum((y_pred-mean(y)).^2);
        ss_error=sum((y-y_pred).^2);
        frat=(ss_model/dim)/(ss_error/(n-dim-1));
        p=1-fcdf(frat,dim,n-dim-1);
        Rsquared=corr(y,y_pred).^2;
        disp(sprintf('   p=%6.4f, F=%8.4f, R^2=%6.4f, recalc',p,frat,Rsquared));
        %
        y_pred_drop=zeros(n,1);
        for k=1:n
            i_drop=setdiff([1:n],n);
            x_drop=x(i_drop,:);
            y_drop=y(i_drop,:); 
            b_drop=regress(y_drop,[ones(n-1,1),x_drop]);
            y_pred_drop(k)=[1 x(k,:)]*b_drop;
        end
        %to do: compare residual error with intrinsic error in behavior;
        %graphics
    end
end
