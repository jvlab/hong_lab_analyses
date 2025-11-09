%hlid_rastim_mds_coords_subsamp: read calcium imaging data from Hong Lab, 
% examine statistics of cross-prep consistency between several transformations of the z-scores
% by doing subsamples of each set of preps
%knit_stats
% This runs on the quantity r saved by hlid_rastim_mds_coords_demo
%
%  See also:  HLID_LOCALOPTS, HLID_RASTIM_COORDS_SUBSAMP_DEMO,
% HLID_RASTIM_MDS_COORDS_SUMM,HLID_RASTIM_MDS_COORDS_SUMM.
%
hlid_opts=hlid_localopts; %set up read_opts and plot_opts 
%
nfiles=length(r.filenames_short);
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

if ~exist('neigs_show_max') neigs_show_max=10; end
%
% for imeth=1:nmeths
%     for submean=0:if_submean
%         meth_text=cat(2,meths{imeth}.name_short,sprintf(' submean=%1.0f',submean))
%         disp(' ');
%         disp(sprintf('processing %s',meth_text));
%         disp(' ');
%         data_in=r.data{imeth,1+submean};
%         [data_align,aux_align]=rs_align_coordsets(data_in,aux);
%         [data_knit,aux_knit]=rs_knit_coordsets(data_align,aux);
%         if (if_write)
%             data_knit.sets{1}.pipeline.embedding_method=setfield(meths{imeth},'if_submean',submean);
%             rs_write_coorddata(filenames_out{imeth,1+submean},data_knit);
%         end
%         r.data_knit{imeth,1+submean}=data_knit;
%         r.knit_stats{imeth,1+submean}=aux_knit.knit_stats;
%         %do stats if nshuffs>0
%         if (nshuffs>0)
%             knit_stats{imeth,1+submean}=aux_knit.knit_stats;
%             knit_stats_setup=aux_knit.knit_stats_setup;
%             knit_stats_setup.dataset_labels=strrep(filenames_short,'.mat','');
%             knit_stats_setup.stimulus_labels=stim_labels;
%             psg_knit_stats_plot(knit_stats{imeth,1+submean},knit_stats_setup);
%             axes('Position',[0.5,0.07,0.01,0.01]); %for text
%             text(0,0,meth_text,'Interpreter','none','FontSize',8);
%             axis off;
%             set(gcf,'Name',sprintf('consensus stats: %s',meth_text));
%         end
%     end %submean
% end %imeth
