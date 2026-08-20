%hlid_vi_pcafilt_summ: summary of KC volumetric imaging, data from George Barnum, Hong Lab,
%
% makes use of output of from hlid_vi_pca_filt_auto.
% 
%   See also:  HLID_VI_PCA_FILT_AUTO, FDR.
%
rmfields={'corrs_avg','corrs_indiv'}; %un-needed large fields to remove
%
if ~exist('pcrits') pcrits=[0.2 0.1 0.05 0.025 0.01]; end %critical values for individual pc's
if ~exist('psumcrits') psumcrits=[1 0.5]; end %critical value for sum of probabilities
%
if ~exist('data_path') data_path='.\'; end
if ~exist('result_file') result_file='hlid_vi_pcafilt_auto_20Aug26.mat'; end
result_file=getinp('result file name','s',[],result_file);
load(cat(2,data_path,filesep,result_file),'results');
%
n_files=length(results);
for ifile=1:n_files
    disp(sprintf('results file entry %2.0f is from %s',ifile,results{ifile}.data_file));
    for fn=1:length(rmfields)
        if isfield(results{ifile},rmfields{fn})
            results{ifile}=rmfield(results{ifile},rmfields{fn});
        end
    end %fn
end
%
strats=struct;
strats.names=cell(0);
strats.type=cell(0);
strats.val=[];
for icrit=1:length(pcrits)
    strats.names{end+1}=sprintf('raw p(F)<%5.3f',pcrits(icrit));
    strats.type{end+1}='raw p';
    strats.val{end+1}=pcrits(icrit);
    strats.names{end+1}=sprintf('fdr p(F)<%5.3f',pcrits(icrit));
    strats.type{end+1}='fdr p';
    strats.val{end+1}=pcrits(icrit);
end
for icrit=1:length(psumcrits)
    strats.names{end+1}=sprintf('max pcs, sum(p)<%5.3f',psumcrits(icrit));
    strats.type{end+1}='max pcs';
    strats.val{end+1}=psumcrits(icrit);
    strats.names{end+1}=sprintf('max pwr, sum(p)<%5.3f',psumcrits(icrit));
    strats.type{end+1}='max pwr';
    strats.val{end+1}=1;
end
n_strats=length(strats.names)
%
for k=1:n_strats
    disp(sprintf(' strategy %3.0f: %s',k,strats.names{k}));
end
%
%process
%
used=cell(1,n_files);
for ifile=1:n_files
    frats_pvals=results{ifile}.frats_pvals;
    eiv_squared=results{ifile}.eiv_squared;
    n_pcs=length(frats_pvals);
    used{ifile}=zeros(n_strats+1,n_pcs);
    for istrat=1:n_strats
        keep=zeros(1,n_pcs);
        switch strats.type{istrat}
            case 'raw p'
                keep=double(frats_pvals<strats.val{istrat});
            case 'fdr p'
                keep=double(frats_pvals<fdr(frats_pvals,strats.val{istrat}));
            case 'max pcs'
                [p_sorted,inds]=sort(frats_pvals);
                limit=find(cumsum(p_sorted)<strats.val{istrat});
                keep(inds(limit))=1;
            case 'max pwr'
        end
        used{ifile}(istrat,:)=keep;
    end %istrat
    used{ifile}(istrat+1,:)=mean(used{ifile}(1:n_strats,:));
    figure;
    set(gcf,'Position',[50 50 800 800]);
    set(gcf,'Name',cat(2,'summary: ',results{ifile}.data_file));
    set(gcf,'NumberTitle','off');
    %
    imagesc(used{ifile}',[0 1]);
    title(result_file,'Interpreter','none');
    ylabel('pc');
    xlabel('strategy');
    set(gca,'XTick',[1:n_strats+1]);
    set(gca,'XTickLabel',[strats.names';'mean']);
    %
    colormap hot;
    colorbar;
    hold on;
    for istrat=0:n_strats
        plot(repmat(0.5+istrat,1,2),[0 n_pcs],'g','LineWidth',2);
    end
    %
    axes('Position',[0.01,0.01,0.01,0.01]);
    text(0,0,results{ifile}.data_file,'Interpreter','none');
    axis off
end %ifile
