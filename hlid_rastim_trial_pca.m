%hlid_rastim_trial_pca: read calcium imaging data from Hong Lab, and analyze at single-trial level
% 
% Produces heatmaps of distances and PCA scree plots and component heatmaps
%   from analyses at stimlus-mean and single-trial level
%   'resid' response is the single-trial response minus the stimulus mean
%   'orth' response is the single-trial response minus the best projection onto the stimulus mean
% Produces results structure.
% Does not mesh with psg routines.
% 
% Analysis is based on z-scored response amplitudes from single trials
%
% Datasets analyzed together must have same sequence of stimuli in
% response_amplitude_stim.stim and response_amplitude_stim.mean_peak,
% and same number of repeats per stimulus (nrepts, defaults to 3)
%
% ROIs with NaNs for all stimuli are kept (these are missing stimuli)
% ROIs with NaNs in any other columns are removed (these are missing ROIs)
% 
% For PCA, number of ROIs without NaNs must be at least as large as number of stimuli.
%
% 21Dec24: fixed pca analysis when stimuli are missing
% 22Dec24: add option for signed comparison of stim with resid and orth, and fix labeling
% 22Dec24: modularize reading files
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO, HLID_RASTIM_TRIAL_READ,
%  HLID_RASTIM_TRIAL_PLOT, HLID_RASTIM_TRIAL_READ, HLID_RASTIM_TRIAL_VIS, HLID_RASTIM_TRIAL_PCA2.
%
hlid_setup;  %invoke hlid_localopts; set up opts_read and opts_plot
%
if ~exist('dist_labels') dist_labels={'Euclidean','1-correl'}; end %can replace by a subset to shorten analysis
ndists=length(dist_labels);
if ~exist('sub_labels') sub_labels={'',' (mean sub)'}; end %can replace  by a subset to shorten analysis
nsubs=length(sub_labels);
%
if ~exist('npca_heatmaps') npca_heatmaps=5; end %number of pc's to show as heatmaps
if ~exist('pca_plotcolors') pca_plotcolors={[0 0 0],[1 0 0],[1 0 1],[0 0 1]}; end %colors for showing pcs
if ~exist('pca_logrange') pca_logrange=2; end %log10 plotting range of variance explained by pcs
if ~exist('npca_maxplot') npca_maxplot=20; end %number of pc's to show in line plot
%
if ~exist('nrepts') nrepts=3; end %number of repeats
%
hlid_rastim_trial_read;
%
%get display order
%
for k=1:nstims
    stimulus_names_display{k}=stimulus_names(k,1:-1+min(find(stimulus_names(k,:)==' ')));
end
fns=fieldnames(display_orders);
if_ok=0;
while if_ok==0
    for k=1:length(fns)
        disp(sprintf(' %1.0f -> display order %s',k,fns{k}));
    end
    k=getinp('choice (0 for native)','d',[0 length(fns)]);
    if k>0
        display_order=display_orders.(fns{k});
    else
        display_order=[];
    end
    if isempty(display_order)
        plot_pos_stim=[1:nstims];
        display_order=stimulus_names_display;
    else
        plot_pos_stim=zeros(1,nstims);
        for k=1:nstims
            posit=min(strmatch(stimulus_names_display{k},display_order,'exact'));
            if posit>0
                plot_pos_stim(k)=posit;
            end
        end
    end
    disp(sprintf('%2.0f stimuli in data files won''t be displayed',sum(double(plot_pos_stim==0))));
    if_ok=getinp('1 if ok','d',[0 1]);
end
ptr_stim=find(plot_pos_stim>0);
ptr_trial=find(reshape(repmat(plot_pos_stim,nrepts,1),1,nstims*nrepts));
%to change plot order for trial-by-trial plot
plot_pos_trial=nrepts*repmat(plot_pos_stim,nrepts,1);
for k=1:nrepts
    plot_pos_trial(k,:)=plot_pos_trial(k,:)+k;
end
plot_pos_trial=(reshape(plot_pos_trial,nstims*nrepts,1))'-nrepts;
%
ticks_stim=[1:nstims];
ticks_trial=[1:nstims]*nrepts-(nrepts-1)/2;
%plot_pos_stim: position to plot each stimulus
%plot_pos_trial: position to plot each trial
%
%
% metadata=cell(nfiles,1);
% dsids=cell(nfiles,1);
% resps_mean=cell(nfiles,1); mean responses within stimuli(nstims,nrois_avail)
% trial_ptrs=cell(nfiles,1); array of (nstims,nrepts)
% resps_trial=cell(nfiles,1); trial-by-trial responses, stimuli unscrambled, (nstims,nrepts,nrois_avail)
% trial_sequence=cell(nfiles,1); stimulus sequence, stimuli as strings
% stims_avail=cell(nfiles,1); list of available stimuli in each file, beginning at 1
% rois_avail=cell(nfiles,1); list of roi numbers kept for analysis, beginning at 1
% rois=cell(nfiles,1):  original rois
% nrois_avail(ifile): number of rois available
%
if_signed_resid_orth=getinp('1 for signed comparisons of resids with resids, and orths with orths','d',[0 1],0);
%
%process distances
%analyze with and without mean-subtraction in each ROI (top row, bot row)
%compute standard distance matrix of mean responses (Euclidean or 1-correlation)
%compute matrix of distances between resids or -resids 
% 
%then do pca of means, and pca of resids, and compare in various ways
%
proc_labels={... %processing methods
    'stims vs stims',...
    'trials vs trials',...
    'resids vs [+/-]resids',...
    'stims vs resids',...
    'orths vs [+/-]orths',...
    'stims vs [+/-]orths'};
if (if_signed_resid_orth==1)
    proc_labels=strrep(proc_labels,'[+/-]','');
end
nprocs=length(proc_labels);
pca_labels={'stims','trials','resids','orths'};
npcas=length(pca_labels);
dists=cell(nfiles,ndists,nsubs,nprocs);
pcas=cell(nfiles,nsubs,npcas);
for ifile=1:nfiles
    for isub=1:nsubs
        resp_stim=resps_mean{ifile};
        resp_trial=resps_trial{ifile};
        switch sub_labels{isub}
            case ''
            case ' (mean sub)'
            resp_xstim=mean(resp_stim,1,'omitnan'); %global mean-
            resp_stim=resp_stim-repmat(resp_xstim,[nstims 1]);
            resp_trial=resp_trial-repmat(reshape(resp_xstim,[1 1 nrois_avail(ifile)]),[nstims nrepts 1]);
        end %isub
        resp_resid=resp_trial-repmat(mean(resp_trial,2),[1 nrepts 1]);
        resp_trial_reshape=reshape(permute(resp_trial,[2 1 3]),[nstims*nrepts,nrois_avail(ifile)]); %dim 1 is now stim (slow) and rept(fast)
        resp_resid_reshape=reshape(permute(resp_resid,[2 1 3]),[nstims*nrepts,nrois_avail(ifile)]);  %dim 1 is now stim (slow) and rept(fast)
        %compute resp_orths_reshape: subtract off best-fitting multiple of mean response
        ind_replic=ceil([1:nstims*nrepts]/nrepts);
        resp_stim_normsq=sum(resp_stim.^2,2);
        resp_trial_dots=sum(resp_trial_reshape.*resp_stim(ind_replic,:),2);
        resp_trial_regress=resp_trial_dots./resp_stim_normsq(ind_replic,:); %multiple of resp_stim to subtract
        resp_orths_reshape=resp_trial_reshape-repmat(resp_trial_regress,[1 size(resp_stim,2)]).*resp_stim(ind_replic,:); %component of trial response orthogonal to mean
        for idist=1:ndists
            switch dist_labels{idist}
                case 'Euclidean'
                    dists{ifile,idist,isub,1}=sqrt(cootodsq(resp_stim,resp_stim));
                    dists{ifile,idist,isub,2}=sqrt(cootodsq(resp_trial_reshape,resp_trial_reshape));
                    dists{ifile,idist,isub,4}=sqrt(cootodsq(resp_stim,resp_resid_reshape));
                    if if_signed_resid_orth==1
                        dists{ifile,idist,isub,3}=sqrt(cootodsq(resp_resid_reshape,resp_resid_reshape));
                        dists{ifile,idist,isub,5}=sqrt(cootodsq(resp_orths_reshape,resp_orths_reshape));
                        dists{ifile,idist,isub,6}=sqrt(cootodsq(resp_stim,resp_orths_reshape));
                    else
                        dists{ifile,idist,isub,3}=sqrt(min(cootodsq(resp_resid_reshape,resp_resid_reshape),cootodsq(resp_resid_reshape,-resp_resid_reshape)));
                        dists{ifile,idist,isub,5}=sqrt(min(cootodsq(resp_orths_reshape,resp_orths_reshape),cootodsq(resp_orths_reshape,-resp_orths_reshape)));
                        dists{ifile,idist,isub,6}=sqrt(min(cootodsq(resp_stim,resp_orths_reshape),cootodsq(resp_stim,-resp_orths_reshape)));
                    end
                case '1-correl'
                    dists{ifile,idist,isub,1}=1-corr(resp_stim',resp_stim');
                    dists{ifile,idist,isub,2}=1-corr(resp_trial_reshape',resp_trial_reshape');
                    dists{ifile,idist,isub,4}=1-corr(resp_stim',resp_resid_reshape');
                    if if_signed_resid_orth==1
                        dists{ifile,idist,isub,3}=1-corr(resp_resid_reshape',resp_resid_reshape');
                        dists{ifile,idist,isub,5}=1-corr(resp_orths_reshape',resp_orths_reshape');
                        dists{ifile,idist,isub,6}=1-corr(resp_stim',resp_orths_reshape');
                    else
                        dists{ifile,idist,isub,3}=1-abs(corr(resp_resid_reshape',resp_resid_reshape'));
                        dists{ifile,idist,isub,5}=1-abs(corr(resp_orths_reshape',resp_orths_reshape'));
                        dists{ifile,idist,isub,6}=1-abs(corr(resp_stim',resp_orths_reshape'));
                    end
            end
        end %idist
        %pca calculations
        pcas{ifile,isub,1}.resp=resp_stim;
        pcas{ifile,isub,2}.resp=resp_trial_reshape;
        pcas{ifile,isub,3}.resp=resp_resid_reshape;
        pcas{ifile,isub,4}.resp=resp_orths_reshape;
        for ipca=1:npcas
%            npcs=min(size(pcas{ifile,isub,ipca}.resp));
            nonans=find(all(~isnan(pcas{ifile,isub,ipca}.resp),2));
            npcs=min(length(nonans),size(pcas{ifile,isub,ipca}.resp,2));
            pcas{ifile,isub,ipca}.npcs=npcs;
            [u,s,v]=svd(pcas{ifile,isub,ipca}.resp(nonans,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
            pcas{ifile,isub,ipca}.s=diag(s(1:npcs,1:npcs));
            pcas{ifile,isub,ipca}.u=nan(size(pcas{ifile,isub,ipca}.resp,1),npcs);
            pcas{ifile,isub,ipca}.u(nonans,:)=u(:,1:npcs);
            pcas{ifile,isub,ipca}.v=v(:,1:npcs);
        end
    end %isub: no mean sub or mean sub
end %ifile
%
%plot distances as heatmaps
%
[nrows_dists,ncols_dists]=nicesubp(nprocs,0.7);
for isub=1:nsubs
    for idist=1:ndists
        for ifile=1:nfiles+1
            if ifile<=nfiles
                tstring=sprintf('file %s, dist: %s %s',dsids{ifile},dist_labels{idist},sub_labels{isub});
            else
                tstring=sprintf('average across %1.0f files, %s to %s, dist: %s %s',nfiles,dsids{1},dsids{nfiles},dist_labels{idist},sub_labels{isub});
            end
            figure;
            set(gcf,'Position',[100 100 1400 800]);
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',tstring);
            for iproc=1:nprocs
                subplot(nrows_dists,ncols_dists,iproc);
                if (ifile<=nfiles)
                    dr=dists{ifile,idist,isub,iproc};
                else
                    dra=zeros(size(dists{1,idist,isub,iproc}));
                    for k=1:nfiles
                        dra(:,:,k)=dists{k,idist,isub,iproc};
                    end
                    dr=mean(dra,3,'omitnan');
                end
                hlid_rastim_trial_plot;
                title(cat(2,proc_labels{iproc},sub_labels{isub}));
            end
            axes('Position',[0.01,0.02,0.01,0.01]); %for text
            text(0,0,tstring,'Interpreter','none');
            axis off;
            disp(sprintf('plotted %s',tstring));
        end %ifile
    end %idist
end %isub
%
%plot pca as log variance explained and heatmaps
%
ncols_pca=2;
%left column: log variance, right column: heatmaps, one for each npca
for isub=1:nsubs
    for ifile=1:nfiles
        tstring=sprintf('pca, file %s %s',dsids{ifile},sub_labels{isub});
        figure;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        for if_norm=0:1
            subplot(2,ncols_pca,1+if_norm*ncols_pca);
            for ipca=1:npcas
                s_diag=pcas{ifile,isub,ipca}.s;
                s_plot=min(npca_maxplot,max(find(s_diag/s_diag(1)>=10^(-pca_logrange))));
                if (if_norm)
                    s_diag=s_diag/s_diag(1);
                end
                hp=semilogy([1:s_plot],s_diag(1:s_plot).^2,'k.-');
                hold on;
                set(hp,'Color',pca_plotcolors{1+mod(ipca-1,length(pca_plotcolors))});
            end
            xlabel('pc');
            set(gca,'XLim',[0 npca_maxplot]);
            if (if_norm)
                set(gca,'YLim',[10^(-pca_logrange) 1]);
                ylabel('rel variance');
            else
                set(gca,'YLim',max(get(gca,'YLim'))*[10^(-pca_logrange) 1]);
                ylabel('variance');
            end
                
            legend(pca_labels);
        end
        for ipca=1:npcas
            u=pcas{ifile,isub,ipca}.u;
            for ih=1:npca_heatmaps
                subplot(npcas,2*npca_heatmaps,npca_heatmaps+ih+(ipca-1)*2*npca_heatmaps);
                up=reshape(u(:,ih),size(u,1)/nstims,nstims);
                up_reorder=zeros(size(up)); 
                up_reorder(:,plot_pos_stim(ptr_stim),:)=up(:,ptr_stim);
                imagesc(up_reorder');
                if size(up,1)>1
                    set(gca,'XTick',[1:size(up,1)]);
                    xlabel('repeat')
                else
                    set(gca,'XTick',[]);
                end
                set(gca,'YTick',[1:nstims]);
                set(gca,'YTickLabel',display_order);
                title(sprintf('%s. pc %2.0f',pca_labels{ipca},ih));
            end
        end
        axes('Position',[0.01,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
        disp(sprintf('plotted %s',tstring));
    end %ifile
end %isub
%
%save results
%
results=struct;
results.if_signed_resid_orth=if_signed_resid_orth;
results.nrepts=nrepts;
results.nstims=nstims;
results.nfiles=nfiles;
results.dists_dims={'d1: file, d2: dist label, d3: sub type, d4: proc, then [stim  or stim, rept] x [stim or stim, rept]'};
results.dists=dists;
results.ndists=ndists;
results.dist_labels=dist_labels;
results.nprocs=nprocs;
results.proc_labels=proc_labels;
results.pcas_dims={'d1: file, d2: sub type, d3: pca type; resp=u*s*v'};
results.pcas=pcas;
results.npcas=npcas;
results.pca_labels=pca_labels;
results.nsubs=nsubs;
results.sub_labels=sub_labels;
results.metadata=metadata;
results.dsids=dsids;
results.stims_avail=stims_avail; %list of available stimuli in each file, beginning at 1
results.rois_avail=rois_avail;
results.rois=rois;
%
results.display_order=display_order;
results.plot_pos_stim=plot_pos_stim;
results.ptr_stim=ptr_stim;
results.ticks_stim=ticks_stim;
results.plot_pos_trial=plot_pos_trial;
results.ptr_trial=ptr_trial;
results.ticks_trial=ticks_trial;
%
disp('results structure created.');
