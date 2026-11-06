% hlid_vi_coords_behav_auto_rs: read volumetric imaging coordinates set files, compare with behavioral data from George Barnum
%
% to do: add r2 and stats to figure; automatic figure saving with dimension, coord and beh data sources
%
% based on hlid_vi_coords_behav_rs, but
%  * automates reading of multiple dimensionally-reduced behavioral indices from George Barnum
%  * directly reads csv files
%
% behavioral data (valence) indicated by color
% reliability of behavior across flies indicated by size (large: less variable)
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, RS_EXTRACT_COORDSETS, HLID_VI_COORDS_KNIT_RS, HLID_VI_COORDS_KNIT_RS_AUTO, REGRESS,
%  HLID_VI_COORDS_BEHAV_RS, HLID_VI_COORDS_BEHAV_DO.
%
hlid_setup;
%
if ~exist('path_beh') path_beh='C:\Users\jdvicto\OneDrive - Weill Cornell Medicine\CloudStorage\From_HongLab\HongLabOrig_for_jdv\volumetric_KC\OdorBehavioralComponents'; end
if ~exist('subdirs_beh') subdirs_beh={'pca_fly_trial_mean','pca_trial_mean_fly_centered','tca_rank_2\fly_odor_weights','tca_rank_3\fly_odor_weights'}; end
if ~exist('fly_csv_string') fly_csv_string='_per_fly.csv';end
%
nsubdirs=length(subdirs_beh);
nbehfiles=0;
nskipped=0;
behfiles_short=cell(0);
data_behs=[];
for isub=1:nsubdirs
    disp(sprintf('behavioral indices in folder %s:',subdirs_beh{isub}));
    path_beh_subdir=cat(2,path_beh,filesep,subdirs_beh{isub},filesep);
    dir_cmd=cat(2,'dir "',path_beh_subdir,filesep,'*',fly_csv_string,'" /b');
    [dos_err,dir_out]=dos(dir_cmd);
    disp(dir_out);
    if dos_err==0
        while ~isempty(dir_out)
            [beh_file,dir_out]=strtok(dir_out);
            if ~isempty(beh_file)
                warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
                table_beh=readtable(cat(2,path_beh_subdir,filesep,beh_file));
                warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
                %verify contents
                if (nbehfiles==0)
                    typenames_beh=table_beh.Properties.VariableDescriptions(2:end)'; %skip the index column
                    nstims_beh=length(typenames_beh);
                    rownames_beh=table_beh.fly;
                    rownums_beh=strmatch('202',rownames_beh); %fly names are 2024_*, etc
                end
                typenames_beh_check=table_beh.Properties.VariableDescriptions(2:end)'; %skip the index column
                nstims_beh_check=length(typenames_beh);
                rownames_beh_check=table_beh.fly;
                rownums_beh_check=strmatch('202',rownames_beh_check); %fly names are 2024_*, etc
                if_ok=1;
                    if nstims_beh_check~=nstims_beh
                        if_ok=0;
                    elseif any(strcmp(typenames_beh,typenames_beh_check)==0) 
                        if_ok=0;
                    end
                    if nbehfiles==5
                        nstims_beh_check=12;
                    end
                    if rownums_beh_check~=rownums_beh
                        if_ok=0;
                    elseif any(strcmp(rownames_beh_check,rownames_beh)==0)
                        if_ok=0;
                    end
                    if if_ok==0
                        disp(sprintf(' mismatch of stims or preps in behavioral file %s; skipped',beh_file));
                        nskipped=nskipped+1;
                    else
                        nbehfiles=nbehfiles+1;
                        behfiles_short{1,nbehfiles}=cat(2,subdirs_beh{isub},' ',beh_file);
                        data_behs(:,:,nbehfiles)=table2array(table_beh(rownums_beh,2:end))'; %rows are stimuli, cols are preps, d3 is behavior index
                    end
             end
        end
    end %valid csv files are present
end %subdirectory
disp(sprintf(' %3.0f valid files of behavioral data files read; %3.0f skipped',nbehfiles,nskipped));
behfiles_short=strrep(behfiles_short,fly_csv_string,'');
behfiles_short=strrep(behfiles_short,'pca_','');
behfiles_short=strrep(behfiles_short,'PC','pc');
behfiles_short=strrep(behfiles_short,'tca_rank_','tca');
behfiles_short=strrep(behfiles_short,' comp',' tc');
behfiles_short=strrep(strrep(behfiles_short,'/','-'),'\','-');
behfiles_short=strrep(behfiles_short,' ','-');
behfiles_short=strrep(behfiles_short,'centered','cent');
behfiles_short=strrep(behfiles_short,'weights','wts');
for ibeh=1:nbehfiles
    disp(sprintf(' index %3.0f -> %s',ibeh,behfiles_short{ibeh}));
end
nstims_beh=size(data_behs,1);
npreps_beh=size(data_behs,2);
disp(sprintf('behavioral data read for %2.0f stimuli from %2.0f preps',nstims_beh,npreps_beh));
%
if ~exist('tol') tol=10^-5; end %for matching stat calcs with matlab
if ~exist('dmax') dmax=5; end
if ~exist('dims_plot') dims_plot=[2 3];end
if ~exist('line_width') line_width=2; end
if ~exist('marker_size') marker_size=12; end
if ~exist('marker_size_nobeh') marker_size_nobeh=8; end
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
disp(sprintf(' %3.0f coordinate sets read.',nsets));
for ibeh=1:nbehfiles
    disp('*******');
    disp(sprintf(' analyzing index %3.0f -> %s',ibeh,behfiles_short{ibeh}));
    disp(' ');
    data_beh=data_behs(:,:,ibeh);
    beh_string=behfiles_short{ibeh};
    mean_beh=mean(data_beh,2);
    %
    %calculate and plot
    hlid_vi_coords_behav_do;
    %
end %ibeh