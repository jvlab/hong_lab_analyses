%hlid_rastim_trial_vis2: auxiliary analyses of single-trial level data
%  
% Run this after hlid_rastim_trial_vis, or after loading results structure.
% See hlid_rastim_trial_vis for documentation and details.
%
%  See also:  HLID_RASTIM_TRIAL_VIS, HLID_RASTIM_TRIAL_PLOT.
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

% for isub=1:nsubs
%     for idist=1:ndists
%         for ifile=nfiles+1:nfiles+1 %only show averages across files
%             if ifile<=nfiles
%                 tstring=sprintf('file %s, dist: %s %s',dsids{ifile},dist_labels{idist},sub_labels{isub});
%             else
%                 tstring=sprintf('average across %2.0f files, %s to %s, dist: %s %s',nfiles,dsids{1},dsids{nfiles},dist_labels{idist},sub_labels{isub});
%             end
%             figure;
%             set(gcf,'Position',[100 100 1400 800]);
%             set(gcf,'NumberTitle','off');
%             set(gcf,'Name',tstring);
%             %
%             axes('Position',[0.01,0.02,0.01,0.01]); %for text
%             text(0,0,tstring,'Interpreter','none');
%             axis off;
%         end %ifile
%     end %idist
% end %isub
