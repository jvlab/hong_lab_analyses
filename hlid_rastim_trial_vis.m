%hlid_rastim_trial_pca: visualize single-trial level data in representational space
%
%  See also:  HLID_SETUP, HLID_LOCALOPTS, HLID_RASTIM2COORDS_DEMO,
%  HLID_RASTIM_TRIAL_READ, HLID_RASTIM_TRIAL_PCA.
%
hlid_setup;  %invoke hlid_localopts; set up read_opts and plot_opts
%
if ~exist('sub_labels') sub_labels={'',' (mean sub)'}; end %can replace  by a subset to shorten analysis
if ~exist('dist_labels') dist_labels={'Euclidean','1-correl'}; end %can replace by a subset to shorten analysis
ndists=length(dist_labels);
%
nsubs=length(sub_labels);
%
if ~exist('nrepts') nrepts=3; end %number of repeats
%
hlid_rastim_trial_read;
%
% metadata=cell(nfiles,1);
% dsids=cell(nfiles,1);
% resps_mean=cell(nfiles,1); %mean responses across stimuli
% trial_ptrs=cell(nfiles,1); %array of (nstims,nrepts)
% resps_trial=cell(nfiles,1); %trial-by-trial responses, stimuli unscrambled, (nstims,nrepts,nrois_avail)
% trial_sequence=cell(nfiles,1); %stimulus sequence, stimuli as strings
% stims_avail=cell(nfiles,1); %list of available stimuli in each file, beginning at 1
% rois_avail=cell(nfiles,1); %list of roi numbers kept for analysis, beginning at 1
% rois=cell(nfiles,1): % original rois

% ma also want to do MDS on the cosine distances
%
%get dimensionality
%get display order

% opts = 
%   struct with fields:
% 
%                  color_norays: [0.2000 0.3000 0.4000]
%                    label_sets: 1
%                  connect_list: [3×2 double]
%     color_norays_connect_mode: 2
%     color_connect_sets_norays: {'r'  'g'  'b'}
%                   if_use_rays: 0
%                    label_list: 'typenames'
% >> save junk
% >> ou=psg_plotcoords(coords,dim_select,sa,rays,opts);
% >> ou=psg_plotcoords(coords,dim_select,sa,rays,opts);
% 131 opts.plot_range=[];
% K>> opts.label_list
% ans =
%   1×6 cell array
%     {'a'}    {'b'}    {'c'}    {'d'}    {'e'}    {'f'}
% K>> opts.label_list
% ans =
%   1×6 cell array
%     {'a'}    {'b'}    {'c'}    {'d'}    {'e'}    {'f'}
% K>> ipt
% ipt =
%      2
% >> save junk
% >> opts.connect_list
% ans =
%      1     2
%      1     3
%      1     4
% >> sa
% sa = 
%   struct with fields:
% 
%     typenames: {'a'  'b'  'c'  'd'  'e'  'f'}
% >> rays
% rays = 
%   struct with no fields.
% >> dim_select
% dim_select =
%      1     2     3
% >>

% %
% ticks_stim=[1:nstims];
% ticks_trial=[1:nstims]*nrepts-(nrepts-1)/2;
% %plot_pos_stim: position to plot each stimulus
% %plot_pos_trial: positoin to plot each trial
% % metadata=cell(nfiles,1);
% % das=cell(nfiles,1);
% % dsids=cell(nfiles,1);
% % resps_mean=cell(nfiles,1); %mean responses across stimuli
% % trial_ptrs=cell(nfiles,1); %array of (nstims,nrepts)
% % resps_trial=cell(nfiles,1); %trial-by-trial responses, stimuli unscrambled, (nstims,nrepts,nrois_avail)
% % trial_sequence=cell(nfiles,1); %stimulus sequence, stimuli as strings
% % stims_avail=cell(nfiles,1); %list of available stimuli in each file, beginning at 1
% % rois_avail=cell(nfiles,1); %list of roi numbers kept for analysis, beginning at 1
% % rois=cell(nfiles,1): % original rois
% %
% %sort trials according to stimuli and summarize
% for ifile=1:nfiles
%     for istim=1:nstims
%         resps_trial{ifile}(istim,:,:)=reshape(das{ifile}.response_amplitude_trials.mean_peak(trial_ptrs{ifile}(istim,:),rois_avail{ifile}),[1 nrepts nrois_avail(ifile)]);
%     end
%     check_mean=max(max(abs(resps_mean{ifile}-squeeze(mean(resps_trial{ifile},2)))));
%     disp(sprintf(' file %2.0f (%30s):  mean responses: %4.0f x %4.0f, stims available: %4.0f, compare computed trial mean and supplied mean: %9.7f',...
%         ifile,dsids{ifile},size(resps_mean{ifile}),length(stims_avail{ifile}),check_mean));
% end
% clear das
% %
% if_signed_resid_orth=getinp('1 for signed comparisons of resids with resids, and orths with orths','d',[0 1],0);
% %
% %process distances
% %analyze with and without mean-subtraction in each ROI (top row, bot row)
% %compute standard distance matrix of mean responses (Euclidean or 1-correlation)
% %compute matrix of distances between resids or -resids 
% % 
% %then do pca of means, and pca of resids, and compare in various ways
% %
% proc_labels={... %processing methods
%     'stims vs stims',...
%     'trials vs trials',...
%     'resids vs [+/-]resids',...
%     'stims vs resids',...
%     'orths vs [+/-]orths',...
%     'stims vs [+/-]orths'};
% if (if_signed_resid_orth==1)
%     proc_labels=strrep(proc_labels,'[+/-]','');
% end
% nprocs=length(proc_labels);
% pca_labels={'stims','trials','resids','orths'};
% npcas=length(pca_labels);
% dists=cell(nfiles,ndists,nsubs,nprocs);
% pcas=cell(nfiles,nsubs,npcas);
% for ifile=1:nfiles
%     for isub=1:nsubs
%         resp_stim=resps_mean{ifile};
%         resp_trial=resps_trial{ifile};
%         if (isub==2)
%             resp_xstim=mean(resp_stim,1,'omitnan'); %global mean-
%             resp_stim=resp_stim-repmat(resp_xstim,[nstims 1]);
%             resp_trial=resp_trial-repmat(reshape(resp_xstim,[1 1 nrois_avail(ifile)]),[nstims nrepts 1]);
%         end %isub
%         resp_resid=resp_trial-repmat(mean(resp_trial,2),[1 nrepts 1]);
%         resp_trial_reshape=reshape(permute(resp_trial,[2 1 3]),[nstims*nrepts,nrois_avail(ifile)]); %dim 1 is now stim (slow) and rept(fast)
%         resp_resid_reshape=reshape(permute(resp_resid,[2 1 3]),[nstims*nrepts,nrois_avail(ifile)]);  %dim 1 is now stim (slow) and rept(fast)
%         %compute resp_orths_reshape: subtract off best-fitting multiple of mean response
%         ind_replic=ceil([1:nstims*nrepts]/nrepts);
%         resp_stim_normsq=sum(resp_stim.^2,2);
%         resp_trial_dots=sum(resp_trial_reshape.*resp_stim(ind_replic,:),2);
%         resp_trial_regress=resp_trial_dots./resp_stim_normsq(ind_replic,:); %multiple of resp_stim to subtract
%         resp_orths_reshape=resp_trial_reshape-repmat(resp_trial_regress,[1 size(resp_stim,2)]).*resp_stim(ind_replic,:); %component of trial response orthogonal to mean
%         for idist=1:ndists
%             switch dist_labels{idist}
%                 case 'Euclidean'
%                     dists{ifile,idist,isub,1}=sqrt(cootodsq(resp_stim,resp_stim));
%                     dists{ifile,idist,isub,2}=sqrt(cootodsq(resp_trial_reshape,resp_trial_reshape));
%                     dists{ifile,idist,isub,4}=sqrt(cootodsq(resp_stim,resp_resid_reshape));
%                     if if_signed_resid_orth==1
%                         dists{ifile,idist,isub,3}=sqrt(cootodsq(resp_resid_reshape,resp_resid_reshape));
%                         dists{ifile,idist,isub,5}=sqrt(cootodsq(resp_orths_reshape,resp_orths_reshape));
%                         dists{ifile,idist,isub,6}=sqrt(cootodsq(resp_stim,resp_orths_reshape));
%                     else
%                         dists{ifile,idist,isub,3}=sqrt(min(cootodsq(resp_resid_reshape,resp_resid_reshape),cootodsq(resp_resid_reshape,-resp_resid_reshape)));
%                         dists{ifile,idist,isub,5}=sqrt(min(cootodsq(resp_orths_reshape,resp_orths_reshape),cootodsq(resp_orths_reshape,-resp_orths_reshape)));
%                         dists{ifile,idist,isub,6}=sqrt(min(cootodsq(resp_stim,resp_orths_reshape),cootodsq(resp_stim,-resp_orths_reshape)));
%                     end
%                 case '1-correl'
%                     dists{ifile,idist,isub,1}=1-corr(resp_stim',resp_stim');
%                     dists{ifile,idist,isub,2}=1-corr(resp_trial_reshape',resp_trial_reshape');
%                     dists{ifile,idist,isub,4}=1-corr(resp_stim',resp_resid_reshape');
%                     if if_signed_resid_orth==1
%                         dists{ifile,idist,isub,3}=1-corr(resp_resid_reshape',resp_resid_reshape');
%                         dists{ifile,idist,isub,5}=1-corr(resp_orths_reshape',resp_orths_reshape');
%                         dists{ifile,idist,isub,6}=1-corr(resp_stim',resp_orths_reshape');
%                     else
%                         dists{ifile,idist,isub,3}=1-abs(corr(resp_resid_reshape',resp_resid_reshape'));
%                         dists{ifile,idist,isub,5}=1-abs(corr(resp_orths_reshape',resp_orths_reshape'));
%                         dists{ifile,idist,isub,6}=1-abs(corr(resp_stim',resp_orths_reshape'));
%                     end
%             end
%         end %idist
%         %pca calculations
%         pcas{ifile,isub,1}.resp=resp_stim;
%         pcas{ifile,isub,2}.resp=resp_trial_reshape;
%         pcas{ifile,isub,3}.resp=resp_resid_reshape;
%         pcas{ifile,isub,4}.resp=resp_orths_reshape;
%         for ipca=1:npcas
% %            npcs=min(size(pcas{ifile,isub,ipca}.resp));
%             nonans=find(all(~isnan(pcas{ifile,isub,ipca}.resp),2));
%             npcs=min(length(nonans),size(pcas{ifile,isub,ipca}.resp,2));
%             pcas{ifile,isub,ipca}.npcs=npcs;
%             [u,s,v]=svd(pcas{ifile,isub,ipca}.resp(nonans,:)); %resp=u*s*v', with u and v both orthogonal, so u*s=resp*v
%             pcas{ifile,isub,ipca}.s=diag(s(1:npcs,1:npcs));
%             pcas{ifile,isub,ipca}.u=nan(size(pcas{ifile,isub,ipca}.resp,1),npcs);
%             pcas{ifile,isub,ipca}.u(nonans,:)=u(:,1:npcs);
%             pcas{ifile,isub,ipca}.v=v(:,1:npcs);
%         end
%     end %isub: no mean sub or mean sub
% end %ifile
% %
% %plot distances as heatmaps
% %
% [nrows_dists,ncols_dists]=nicesubp(nprocs,0.7);
% for isub=1:nsubs
%     for idist=1:ndists
%         for ifile=1:nfiles+1
%             if ifile<=nfiles
%                 tstring=sprintf('file %s, dist: %s %s',dsids{ifile},dist_labels{idist},sub_labels{isub});
%             else
%                 tstring=sprintf('average across files %s to %s, dist: %s %s',dsids{1},dsids{nfiles},dist_labels{idist},sub_labels{isub});
%             end
%             figure;
%             set(gcf,'Position',[100 100 1400 800]);
%             set(gcf,'NumberTitle','off');
%             set(gcf,'Name',tstring);
%             for iproc=1:nprocs
%                 subplot(nrows_dists,ncols_dists,iproc);
%                 if (ifile<=nfiles)
%                     dr=dists{ifile,idist,isub,iproc};
%                 else
%                     dra=zeros(size(dists{1,idist,isub,iproc}));
%                     for k=1:nfiles
%                         dra(:,:,k)=dists{k,idist,isub,iproc};
%                     end
%                     dr=mean(dra,3,'omitnan');
%                 end
%                 hlid_rastim_trial_plot;
%                 title(cat(2,proc_labels{iproc},sub_labels{isub}));
%             end
%             axes('Position',[0.01,0.02,0.01,0.01]); %for text
%             text(0,0,tstring,'Interpreter','none');
%             axis off;
%             disp(sprintf('plotted %s',tstring));
%         end %ifile
%     end %idist
% end %isub
% %
% %plot pca as log variance explained and heatmaps
% %
% ncols_pca=2;
% %left column: log variance, right column: heatmaps, one for each npca
% for isub=1:nsubs
%     for ifile=1:nfiles
%         tstring=sprintf('pca, file %s %s',dsids{ifile},sub_labels{isub});
%         figure;
%         set(gcf,'Position',[100 100 1200 800]);
%         set(gcf,'NumberTitle','off');
%         set(gcf,'Name',tstring);
%         for if_norm=0:1
%             subplot(2,ncols_pca,1+if_norm*ncols_pca);
%             for ipca=1:npcas
%                 s_diag=pcas{ifile,isub,ipca}.s;
%                 s_plot=min(npca_maxplot,max(find(s_diag/s_diag(1)>=10^(-pca_logrange))));
%                 if (if_norm)
%                     s_diag=s_diag/s_diag(1);
%                 end
%                 hp=semilogy([1:s_plot],s_diag(1:s_plot).^2,'k.-');
%                 hold on;
%                 set(hp,'Color',pca_plotcolors{1+mod(ipca-1,length(pca_plotcolors))});
%             end
%             xlabel('pc');
%             set(gca,'XLim',[0 npca_maxplot]);
%             if (if_norm)
%                 set(gca,'YLim',[10^(-pca_logrange) 1]);
%                 ylabel('rel variance');
%             else
%                 set(gca,'YLim',max(get(gca,'YLim'))*[10^(-pca_logrange) 1]);
%                 ylabel('variance');
%             end
% 
%             legend(pca_labels);
%         end
%         for ipca=1:npcas
%             u=pcas{ifile,isub,ipca}.u;
%             for ih=1:npca_heatmaps
%                 subplot(npcas,2*npca_heatmaps,npca_heatmaps+ih+(ipca-1)*2*npca_heatmaps);
%                 up=reshape(u(:,ih),size(u,1)/nstims,nstims);
%                 up_reorder=zeros(size(up)); 
%                 up_reorder(:,plot_pos_stim(ptr_stim),:)=up(:,ptr_stim);
%                 imagesc(up_reorder');
%                 if size(up,1)>1
%                     set(gca,'XTick',[1:size(up,1)]);
%                     xlabel('repeat')
%                 else
%                     set(gca,'XTick',[]);
%                 end
%                 set(gca,'YTick',[1:nstims]);
%                 set(gca,'YTickLabel',display_order);
%                 title(sprintf('%s. pc %2.0f',pca_labels{ipca},ih));
%             end
%         end
%         axes('Position',[0.01,0.02,0.01,0.01]); %for text
%         text(0,0,tstring,'Interpreter','none');
%         axis off;
%         disp(sprintf('plotted %s',tstring));
%     end %ifile
% end %isub
% %
% %save results
% %
% results=struct;
% results.if_signed_resid_orth=if_signed_resid_orth;
% results.nrepts=nrepts;
% results.nstims=nstims;
% results.nfiles=nfiles;
% results.dists_dims={'d1: file, d2: dist label, d3: sub type, d4: proc, then [stim  or stim, rept] x [stim or stim, rept]'};
% results.dists=dists;
% results.ndists=ndists;
% results.dist_labels=dist_labels;
% results.nprocs=nprocs;
% results.proc_labels=proc_labels;
% results.pcas_dims={'d1: file, d2: sub type, d3: pca type; resp=u*s*v'};
% results.pcas=pcas;
% results.npcas=npcas;
% results.pca_labels=pca_labels;
% results.nsubs=nsubs;
% results.sub_labels=sub_labels;
% results.metadata=metadata;
% results.dsids=dsids;
% results.stims_avail=stims_avail; %list of available stimuli in each file, beginning at 1
% results.rois_avail=rois_avail;
% results.rois=rois;
% %
% disp('results structure created.');
