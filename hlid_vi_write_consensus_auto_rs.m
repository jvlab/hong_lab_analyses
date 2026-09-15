% hlid_vi_write_consensus_auto_rs: read volumetric imaging coordinates set files,
% automated knit into consensus and write the consensus files, based on outputs of hlid_vi_pcafilt_coords_auto.
%
% Uses results from hlid_vi_pcafilt_coords_auto to determine names of coordinate files and the options usef for conversion to coordinates
% Reads these coordinate files, aligns, does consensus, and writes
%
%   See also:  HLID_SETUP, RS_GET_COORDSETS, HLID_VI_PCAFILT_COORDS_AUTO,
%   HLID_VI_COORDS_KNIT_RS, HLIS_METHS_DEFINE, RS_GET_COORDSETS, RS_KNIT_COORDSETS, RS_WRITE_COORDSETS.
%
hlid_setup;
%
if ~exist('prefix_remove') prefix_remove='hlid_vi_'; end
%
opts_read=struct();
opts_read.input_type=1; %just data
opts_read.if_auto=1; %no confirmation needed
opts_read.type_class_def='hlid';
opts_read.type_coords_def='zeros';
opts_read.type_class_aux=opts_read.type_class_def;
opts_read.if_log=0;
%
if ~exist('opts_align') opts_align=struct(); end %for psg_align_coordsets
opts_align.if_log=0;
if ~exist('opts_check') opts_check=struct(); end
opts_check.if_warn=0;
%
if ~exist('results_file') results_file='hlid_vi_pcafilt_coords_auto_11Sep26.mat'; end
if ~exist('coord_path') coord_path='./data/kc_vi'; end
results_file=getinp('file with results from hlid_vi_pcafilt_coords_auto','s',[],results_file);
coord_path=getinp('path to coordinate files created by hlid_vi_pcafilt_coords_auto','s',[],coord_path);
coord_path=strrep(strrep(coord_path,'/',filesep),'\',filesep);
%
load(results_file,'results');
results_dims=size(results);
nsets=results_dims(1);
nvariants=prod(results_dims(2:end));
results_res=reshape(results,nsets,nvariants);
disp('loaded results, with dimensions')
disp(results_dims)
disp(sprintf('nsets=%3.0f  nvariants=%5.0f',nsets,nvariants));
disp(sprintf('first coord file of first variant: %s',results_res{1,1}.coord_file));
disp(sprintf(' last coord file of first variant: %s',results_res{end,1}.coord_file));
disp(sprintf('first coord file of  last variant: %s',results_res{1,end}.coord_file));
disp(sprintf(' last coord file of  last variant: %s',results_res{end,end}.coord_file));
%
disp(' ');
for k=1:nsets
    disp(sprintf(' coord file %2.0f of %2.0f: %s',k,nsets,results_res{k,1}.coord_file));
end
%
for k=1:size(results,2)
    sfilts(1,k)=results{1,k,1,1,1,1}.sfilt;
end
for k=1:size(results,4)
    pcrits(1,k)=results{1,1,1,k,1,1}.pcrit;
end
disp('spatial filter values');
disp(sfilts);
disp('pcrit values');
disp(pcrits);
%
if ~exist('dim_max_in') dim_max_in=10; end
%
%determine which methods should be knitted with normalized scales
meths=hlid_meths_define;
meth_strings=cell(1,length(meths));
for imeth=1:length(meths)
    meth_strings{imeth}=meths{imeth}.name_file;
end
%
if_allow_scales=contains(meth_strings,'euc');
if_ok=0;
while (if_ok==0)
    for imeth=1:length(meths)
        disp(sprintf(' method %2.0f is %20s knitted with allow_scale set to %1.0f',imeth,meth_strings{imeth},if_allow_scales(imeth)));
    end
    if_ok=getinp('1 if ok','d',[0 1]);
    if ~if_ok
        if_allow_scales=getinp('new values','d',[0 1],if_allow_scales);
    end
end
%
dim_max_in=getinp('maximum dimension to consider','d',[1 24],dim_max_in);
if_remove_pipeline=getinp('1 to remove pipeline field from results_knit consensus data','d',[0 1],1);
replace_string=getinp('string to replace prep ID with','s',[],'consensus');
%
opts_knit=struct();
opts_knit.allow_reflection=1;
opts_knit.allow_offset=1;
opts_knit.allow_scale=0;
opts_knit.if_normscale=1;
opts_knit.if_pca=0; %if_c2p handled later
opts_knit.max_niters=1000;
opts_knit.pcon_init_method=0;
opts_knit.if_stats=0;
opts_knit.if_plot=0;
opts_knit.dim_max_in=dim_max_in;
opts_knit.nshuffs=0;
opts_knit.if_frozen=1;
opts_knit.if_log=0;
%
warnings=[];
for ivariant=1:nvariants
    fullnames=cell(1,nsets);
    fullnames_out=cell(1,nsets);
    for iset=1:nsets
        fullnames{iset}=cat(2,coord_path,filesep,results_res{iset,ivariant}.coord_file);
        %to form output file name, replace _coords_2025-01-02-a_sf by _coords_consensus_sf
        rep_start=strfind(fullnames{iset},'_coords_')+length('_coords_')-2;
        rep_end=strfind(fullnames{iset},'_sf')+1;
        fullnames_out{iset}=cat(2,fullnames{iset}(1:rep_start),'_',replace_string,'_',fullnames{iset}(rep_end:end));
    end
    aux=struct;
    aux.nsets=nsets;
    [data_read,aux_read_out]=rs_get_coordsets(fullnames,setfield(aux,'opts_read',opts_read));
    disp('***********')
    disp(sprintf('variant %5.0f',ivariant));
    disp(sprintf('first file %s',fullnames{1}));
    disp(sprintf(' last file %s',fullnames{nsets}));
    %
    [data_aligned,aux_align_out]=rs_align_coordsets(data_read,setfields(struct(),{'opts_align','opts_check'},{opts_align,opts_check}));
    nstims_all=data_aligned.sets{1}.nstims;
    disp('aligned');
    %
    meth_string=results_res{iset,ivariant}.meth_string;
    imeth=find(contains(meth_strings,meth_string)>0);
    fullname_out=unique(fullnames_out);
    if length(fullname_out)~=1
        wmsg=sprintf('for variant %4.0f, input file names are not consistent',ivariant);
        warning(wmsg);
        warnings=strvcat(warnings,wmsg);
        warnings=strvcat(warnings,fullnames{:});
    end
    if length(imeth)==1
        opts_knit.allow_scale=if_allow_scales(imeth);
        aux.opts_knit=opts_knit;
        [data_consensus,aux_knit_out]=rs_knit_coordsets(data_aligned,setfields(struct(),{'opts_knit','opts_check'},{opts_knit,opts_check}));
        disp(sprintf('knit with dim_max_in=%3.0f, pcon_init_method=%3.0f, allow_scale=%1.0f',opts_knit.dim_max_in,opts_knit.pcon_init_method,opts_knit.allow_scale));
        if if_remove_pipeline
            data_consensus.sets{1}=struct();
        end
        disp(sprintf('output consensus file will be %s',fullnames_out{1}));
        rs_write_coorddata(fullnames_out{1},data_consensus);
    else
        wmsg=sprintf('for variant %4.0f, method %s not recognized',ivariant,meth_string);
        warning(wmsg);
        warnings=strvcat(warnings,wmsg);
    end
end
if ~isempty(warnings)
    disp(warnings);
end
