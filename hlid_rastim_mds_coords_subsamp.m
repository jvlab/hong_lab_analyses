%hlid_rastim_mds_coords_subsamp: read calcium imaging data from Hong Lab, 
% examine statistics of cross-prep consistency between several transformations of the z-scores
% by doing subsamples of each set of preps and extracting the variance fields from knit_stats
% (but no shuffles are done for the subsamples)
% 
% This runs on the quantity r saved by hlid_rastim_mds_coords_demo,
% and creates r_subsamp
%
%
%  See also:  HLID_LOCALOPTS, HLID_RASTIM_COORDS_SUBSAMP_DEMO,
% HLID_RASTIM_MDS_COORDS_SUMM,HLID_RASTIM_MDS_COORDS_SUMM.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
nfiles=length(r.filenames_short);
nstims=r.nstims;
nmeths=length(r.meths);
if ~exist('subsamp_sizes') subsamp_sizes=[2 3 4 6 8]; end
if ~exist('exhaust_max') exhaust_max=500; end
%
disp(sprintf(' analyzing %2.0f files, %2.0f methods',nfiles,nmeths));
for ifile=1:nfiles
    disp(sprintf('file %2.0f->%s',ifile,r.filenames_short{ifile}));
end
for imeth=1:nmeths
    disp(sprintf('method %2.0f->%s',imeth,r.meths{imeth}.name_full));
end
dim_max_in=getinp('max dimension to use','d',[1 min(r.maxdim(:))]);
%
subsamp_sizes=getinp('subsample list','d',[2 Inf],subsamp_sizes);
nsizes=length(subsamp_sizes);
exhaust_max=getinp('max number of subsamples for exhaustive list','d',[0 Inf],exhaust_max);
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if (if_frozen~=0) 
    rng('default');
    if (if_frozen<0)
        rand(1,abs(if_frozen));
    end
else
    rng('shuffle');
end
%
subsamp_lists=cell(1,nsizes);
for isize=1:nsizes
    subsamp_size=subsamp_sizes(isize);
    if nfiles>=subsamp_size
        exhaust=nchoosek(nfiles,subsamp_size);
        if exhaust<=exhaust_max
            subsamp_lists{isize}=nchoosek([1:nfiles],subsamp_size);
        else
            subsamp_lists{isize}=zeros(exhaust_max,subsamp_size);
            for k=1:exhaust_max
                subsamp_lists{isize}(k,:)=randperm(nfiles,subsamp_size);
            end
        end
    else
        subsamp_lists{isize}=zeros(0,subsamp_size);
    end
    disp(sprintf('for subsample size %2.0f, %5.0f subsamples prepared',subsamp_size,size(subsamp_lists{isize},1)));
end
%
%
%align and knit with rs package
%
aux=struct;
aux.opts_align.if_log=0;
aux.opts_knit.if_log=0; %since we will be edoing a lot of subsamples
aux.opts_knit.allow_scale=0;
aux.opts_knit.if_normscale=0;
aux.opts_knit.keep_details=0; 
aux.opts_knit.if_stats=1; %need statistics of variance
aux.opts_knit.nshuffs=0;
aux.opts_knit.if_plot=0; %no plotting
aux.opts_knit.dim_max_in=dim_max_in;
aux.opts_check.if_warn=0; %no need to check the subsamples
aux_warn=aux;
aux_warn.opts_check.if_warn=1;
%
if_submean=1;
%
r_subsamp.nfiles=nfiles;
r_subsamp.nstims=nstims;
r_subsamp.nmeths=nmeths;
r_subsamp.subsamp_sizes=subsamp_sizes;
r_subsamp.nsizes=nsizes;
r_subsamp.exhaust_max=exhaust_max;
r_subsamp.dim_max_in=dim_max_in;
r_subsamp.subsamp_lists=subsamp_lists;
r_subsamp.stats=cell(nmeths,1+if_submean,nsizes);
r_subsamp.stats_dims={'d1: methods, d2: submean, d3: size of subsample'};
for isize=1:nsizes
    subsamp_size=subsamp_sizes(isize);
    nsubsamps=size(subsamp_lists{isize},1);
    disp(' ')
    disp(sprintf(' subsampling %2.0f files (%5.0f subsamples)',subsamp_size,nsubsamps));
    if nsubsamps>0
        for imeth=1:nmeths
            for submean=0:if_submean
                data_in=r.data{imeth,1+submean}; %data from all files
                %check consistency (first pass only)
                if (imeth==1) & (submean==0)
                    [data_align,aux_align]=rs_align_coordsets(data_in,aux_warn);
                    [data_knit,aux_knit]=rs_knit_coordsets(data_align,aux_warn);
                end
                meth_text=cat(2,r.meths{imeth}.name_short,sprintf(' submean=%1.0f',submean));
                disp(sprintf('processing %s',meth_text));
                sv=struct;
                sv.rmsdev_setwise=NaN(dim_max_in,subsamp_size,nsubsamps); %number of files = size of subsample
                sv.rmsdev_stmwise=NaN(dim_max_in,nstims,nsubsamps);
                sv.rmsdev_overall=NaN(dim_max_in,1,nsubsamps);
                sv.rmsavail_setwise=NaN(dim_max_in,subsamp_size,nsubsamps);
                sv.rmsavail_stmwise=NaN(dim_max_in,nstims,nsubsamps);
                sv.rmsavail_overall=NaN(dim_max_in,1,nsubsamps);
                sv_fields=fieldnames(sv);
                for isubsamp=1:nsubsamps
                    select=subsamp_lists{isize}(isubsamp,:);
                    data_in_sub.ds=data_in.ds(select);
                    data_in_sub.sas=data_in.sas(select);
                    data_in_sub.sets=data_in.sets(select);
                    [data_align_sub,aux_align_sub]=rs_align_coordsets(data_in_sub,aux);
                    [data_knit_sub,aux_knit_sub]=rs_knit_coordsets(data_align_sub,aux);
                    %in the subsample, not all stimuli may have been present
                    sub_typenames=data_knit_sub.sas{1}.typenames;
                    typenames=data_knit.sas{1}.typenames;
                    stm_ptrs=zeros(length(sub_typenames),1);
                    for istim=1:length(sub_typenames)
                        stm_ptrs(istim)=strmatch(sub_typenames{istim},typenames,'exact');
                    end
                    for isv=1:length(sv_fields)
                        sv_field=sv_fields{isv};
                        if ~isempty(strfind(sv_field,'_stm')) %is  it a tally by stimuli?
                            sv.(sv_field)(:,stm_ptrs,isubsamp)=aux_knit_sub.knit_stats.(sv_field);
                        else
                            sv.(sv_field)(:,:,isubsamp)=aux_knit_sub.knit_stats.(sv_field);
                        end
                    end
                end
                r_subsamp.stats{imeth,1+submean,isize}=sv;
            end %submean
        end %imeth
    end %nsubsamps>0
end %isize
