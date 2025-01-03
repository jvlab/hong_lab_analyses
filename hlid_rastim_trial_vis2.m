%hlid_rastim_trial_vis2: auxiliary analyses of single-trial level data
%  
% Run this after hlid_rastim_trial_vis, or after loading results structure.
% See hlid_rastim_trial_vis for documentation and details.
%
% degrees of freedom assumes that stimuli are either completely missing from a file, or have
% correct number of repeats
% 
% 02Jan25: add option for recentering residuals prior to covariance calculation
%
%  See also:  HLID_RASTIM_TRIAL_VIS, HLID_RASTIM_TRIAL_PLOT, RANDORTHU, RANDORTHU_GEN.
%
nrepts=results.nrepts;
nstims=results.nstims;
nfiles=results.nfiles;
nsubs=results.nsubs;
sub_labels=results.sub_labels;
npreprocs=results.npreprocs;
preproc_labels=results.preproc_labels;
nsts=results.nsts;
st_labels=results.st_labels;
dmax=results.dmax;
coords_stim=results.coords_stim;
coords_trial=results.coords_trial;
metadata=results.metadata;
dsids=results.dsids;
stims_avail=results.stims_avail; %list of available stimuli in each file, beginning at 1
rois_avail=results.rois_avail;
rois=results.rois;
%
coords_consensus=results.coords_consensus;
xforms_consensus=results.xforms_consensus;
coords_devs_xform=results.coords_devs_xform;
% results.coords_devs_xform_dims={'{d1: 1(consensus), d2: sub type, d3: preproc type, d4: stim or trial}{model dim}(stim,coord,rep,file)'};
%
nrepts_tot=nrepts*nfiles;
disp(sprintf('analyzing covariances from %3.0f responses per stimulus (%2.0f datasets, %2.0f repeats per dataset)',nrepts_tot,nfiles,nrepts))
coveigs=cell(nsubs,npreprocs,nsts);
coveigs_surr=cell(nsubs,npreprocs,nsts);
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if ~exist('ndraws') ndraws=100; end
ndraws=getinp('ndraws','d',[1 10000],ndraws);
dmax_use=getinp('max dim to calculate','d',[1 dmax]);
if_recenter=getinp('1 to recenter resids prior to covariance','d',[0 1],0);
surrtypes={'random direction','keep first eiv, then random'};
ntypes=length(surrtypes); %number of surrogate types
if ~exist('quantiles_show') quantiles_show=[.01 .05 .5 .95 .99]; end
%
for isub=1:nsubs
    for ipreproc=1:npreprocs
        for ist=1:nsts
                disp(' ')
                tstring=sprintf('merging across %2.0f files, %s to %s, %s %s, pc on %s, recenter for covariances: %1.0f',...
                    nfiles,dsids{1},dsids{nfiles},sub_labels{isub},preproc_labels{ipreproc},st_labels{ist},if_recenter);               
                disp(tstring);
                coveigs{isub,ipreproc,ist}=cell(1,dmax_use);
                coveigs_surr{isub,ipreproc,ist}=cell(1,dmax_use);
                for  idim=1:dmax_use
                     if (if_frozen~=0) 
                        rng('default');
                        if (if_frozen<0)
                            rand(1,abs(if_frozen));
                        end
                    else
                        rng('shuffle');
                    end
                    coveigs{isub,ipreproc,ist}{idim}=zeros(nstims,idim);
                    coveigs_surr{isub,ipreproc,ist}{idim}=zeros(nstims,idim,ndraws,ntypes);
                    coords_use=reshape(coords_devs_xform{1,isub,ipreproc,ist}{idim},[nstims,idim,nrepts_tot]);
                    % for each stimulus, use coords_use to get a crude
                    % estimate of the covariances of the residuals
                    for istim=1:nstims
                        coords_cov=reshape(coords_use(istim,:,:),[idim nrepts_tot]); %coordinates to use for covariance estimate
                        coords_cov=coords_cov(:,all(~isnan(coords_cov),1)); %remove NaNs
                        if (if_recenter)
                            coords_cov=coords_cov-repmat(mean(coords_cov,2),[1 size(coords_cov,2)]);
                        end
                        cov_dof=size(coords_cov,2)*(nrepts-1)/nrepts; %lose one degree of freedom for each file in which data (nrepts) are present
                        cov_est=coords_cov*coords_cov'/cov_dof; %denominator is degrees of freedom 
                        %now look at eigenvalue of cov_est, and accumulate across stimuli
                        [eivecs,eivals]=eig(cov_est);
                        eivals=diag(eivals)';
                        [coveigs{isub,ipreproc,ist}{idim}(istim,:),idx]=sort(eivals(:),'descend');
                        eivecs=eivecs(:,idx);
                        %
                        %surrogate type 1: random direction
                        itype=1;
                        coords_cov_surr=zeros(size(coords_cov));
                        for idraw=1:ndraws
                            for irept=1:size(coords_cov,2)
                                rotm=randorthu(idim); %rotate each residual by a different random unitary matrix
                                coords_cov_surr(:,irept)=rotm*coords_cov(:,irept);
                            end
                            cov_est_surr=coords_cov_surr*coords_cov_surr'/cov_dof;
                            eivals_surr=eig(cov_est_surr); 
                            coveigs_surr{isub,ipreproc,ist}{idim}(istim,:,idraw,itype)=sort(eivals_surr(:),'descend')';
                        end                                               %
                        %surrogate type 2: keep first direction, then random 
                        itype=2;
                        coords_cov_surr=zeros(size(coords_cov));
                        for idraw=1:ndraws
                            for irept=1:size(coords_cov,2)
                                rotm=randorthu_gen(idim,eivecs(:,1)); %rotate each residual by a different random unitary matrix that preserves first eigenvector
                                coords_cov_surr(:,irept)=rotm*coords_cov(:,irept);
                                if (if_recenter)
                                    coords_cov_surr=coords_cov_surr-repmat(mean(coords_cov_surr,2),[1 size(coords_cov_surr,2)]);
                                end                              
                            end
                            cov_est_surr=coords_cov_surr*coords_cov_surr'/cov_dof;
                            eivals_surr=eig(cov_est_surr); 
                            coveigs_surr{isub,ipreproc,ist}{idim}(istim,:,idraw,itype)=sort(eivals_surr(:),'descend')';
                        end
                    end %istim
                end %idim
                %
                %plot and display quantiles
                %
                nrows=ntypes;
                ncols=dmax-1;
                figure;
                set(gcf,'Position',[100 100 1400 800]);
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',tstring)
                %
                disp(sprintf('mean of eigenvalues of normalized residuals covariance matrix across %2.0f stimuli',nstims))
                for idim=1:dmax_use
                    ceig=coveigs{isub,ipreproc,ist}{idim};
                    ceig_norm=ceig./repmat(sum(ceig,2),[1 idim]);
                    disp(' ');
                    disp(sprintf(' dim %2.0f:        %s',idim,sprintf(' %7.4f',mean(ceig_norm,1))));
                    %
                    for itype=1:ntypes
                        disp(sprintf(' surrogate %1.0f: %s',itype,surrtypes{itype}));
                        ceig_surr=coveigs_surr{isub,ipreproc,ist}{idim}(:,:,:,itype);
                        ceig_norm_surr=ceig_surr./repmat(sum(ceig_surr,2),[1 idim 1]);
                        disp(sprintf('  means:        %s',sprintf(' %7.4f',mean(mean(ceig_norm_surr,1),3))));
                        quantiles=(sum(repmat(ceig_norm,[1 1 ndraws])>=ceig_norm_surr,3))/ndraws;
                        for iq=1:length(quantiles_show)
                            disp(sprintf(' count q>%5.3f: %s',quantiles_show(iq),sprintf(' %7.0f',sum(quantiles>quantiles_show(iq)))));                     
                        end %iq
                        %
                        if (idim>=2)
                            subplot(nrows,ncols,(idim-1)+(itype-1)*ncols);
                            plot(quantiles','LineWidth',2);
                            set(gca,'XLim',[0.5 idim+0.5]);
                            set(gca,'XTick',[1:1:idim]);
                            xlabel('eig');
                            set(gca,'YLim',[0 1]);
                            ylabel(sprintf('quantile, surr %1.0f',itype));
                            title(sprintf('space dim %1.0f',idim));
                        end
                    end %itype
                end %idim
                %
                axes('Position',[0.01,0.05,0.01,0.01]); %for text
                text(0,0,sprintf('ndraws=%4.0f,  if_frozen=%4.0f',ndraws,if_frozen),'Interpreter','none');
                axis off;
                %
                axes('Position',[0.01,0.02,0.01,0.01]); %for text
                text(0,0,tstring,'Interpreter','none');
                axis off;
        end %ifile
    end %idist
end %isub
