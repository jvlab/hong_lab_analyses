%hlid_rastim_trial_vis2: auxiliary analyses of single-trial level data
%  
% Run this after hlid_rastim_trial_vis, or after loading results structure.
% See hlid_rastim_trial_vis for documentation and details.
%
% To do: add a second surrogate type in which there is a random rotation
% around the response direction; use extorthb to extend a vector to an orthogonal 
% basis, and then generate random rotations with one (or more) axes fixed
%
%  See also:  HLID_RASTIM_TRIAL_VIS, HLID_RASTIM_TRIAL_PLOT, RANDORTHU.
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
nrepts_dof=(nrepts-1)*nfiles;
disp(sprintf('analyzing covariances from %3.0f responses per stimulus (%2.0f d.o.f., %2.0f datasets, %2.0f repeats per dataset)',nrepts_tot,nrepts_dof,nfiles,nrepts))
coveigs=cell(nsubs,npreprocs,nsts);
coveigs_surr=cell(nsubs,npreprocs,nsts);
%
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time, <0 for a specific seed','d',[-10000 1],1);
if ~exist('ndraws') ndraws=100; end
ndraws=getinp('ndraws','d',[1 10000],ndraws);
dmax_use=getinp('max dim to calculate','d',[1 dmax]);
surrtypes={'random direction'};
disp(' surrogate types:');
disp(surrtypes');
ntypes=length(surrtypes); %number of surrogate types
%
for isub=1:nsubs
    for ipreproc=1:npreprocs
        for ist=1:nsts
                tstring=sprintf('average across %2.0f files, %s to %s, %s %s,pc on %s',nfiles,dsids{1},dsids{nfiles},sub_labels{isub},preproc_labels{ipreproc},st_labels{ist});               
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
                        cov_est=coords_cov*coords_cov'/nrepts_dof; %denominator is degrees of freedom 
                        %now look at eigenvalue of cov_est, and accumulate across stimuli
                        eivals=eig(cov_est); 
                        coveigs{isub,ipreproc,ist}{idim}(istim,:)=sort(eivals(:),'descend')';
                        %
                        %surrogate type 1: random direction
                        itype=1;
                        coords_cov_surr=zeros(size(coords_cov));
                        for idraw=1:ndraws
                            for irept=1:nrepts_tot
                                rotm=randorthu(idim); %rotate each residual by a different random unitary matrix
                                coords_cov_surr(:,irept)=rotm*coords_cov(:,irept);
                            end
                            cov_est_surr=coords_cov_surr*coords_cov_surr'/nrepts_dof;
                            eivals_surr=eig(cov_est_surr); 
                            coveigs_surr{isub,ipreproc,ist}{idim}(istim,:,idraw,itype)=sort(eivals_surr(:),'descend')';
                        end
                    end
                end
                %
                disp('mean of eigenvalues of normalized residuals covariance matrix')
                for idim=1:dmax_use
                    ceig=coveigs{isub,ipreproc,ist}{idim};
                    ceig_norm=ceig./repmat(sum(ceig,2),[1 idim]);
                    disp(sprintf(' dim %2.0f: %s',idim,sprintf(' %7.4f',mean(ceig_norm,1))));
                    %
                    for itype=1:ntypes
                        ceig_surr=coveigs_surr{isub,ipreproc,ist}{idim}(:,:,:,itype);
                        ceig_norm_surr=ceig_surr./repmat(sum(ceig_surr,2),[1 idim 1]);
                        disp(sprintf('surr %2.0f: %s',itype,sprintf(' %7.4f',mean(mean(ceig_norm_surr,1),3))));
                    end
                end
                % figure;
                % set(gcf,'Position',[100 100 1400 800]);
                % set(gcf,'NumberTitle','off');
                % set(gcf,'Name',tstring);
                % %
                % axes('Position',[0.01,0.02,0.01,0.01]); %for text
                % text(0,0,tstring,'Interpreter','none');
                % axis off;
        end %ifile
    end %idist
end %isub
