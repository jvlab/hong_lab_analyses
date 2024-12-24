function [opts_used,opts_trial_used]=hlid_rastim_trial_vis_plot(coords_mean,coords_trial,dimlist,colors,labels,opts)
% [opts_used,opts_trial_used]=hlid_rastim_trial_vis_plot(coords_mean,coords_trial,dimlist,colors,labels,opts)
% plots mean responses and signle-trial responses on same axeslsit
%
% coords_mean: mean responses (nstims,dims)
% coords_trial: signle-trial responses (nstims,dims,nrepts)
% dimlist: dimensions to plot
% colors: colors for each trial, as cell array of length nrepts
% labels: stimulus labels, as cell array of length nstims
% opts: options, can be empty
%    axis_handle, if specified, is axis to plot into
%    see fields filled by filldefault for other customizations
%
% opt_used: options used from composite plot
%   opts_used.axis_handle: axis handle
% opts_trial_used: cell array of optoins used for each trial's plot
%
%  See also: HLID_RASTIM_TRIAL_VIS, PSG_PLOTCOORDS, HLID_RASTM_TRIAL_VIS_LEGEND.
%
if nargin<=5
    opts=struct;
end
nrepts=size(coords_trial,3);
opts.if_use_rays=0;
opts.label_sets=1;
opts.if_legend=0; %add legend locally
opts.color_norays_connect_mode=2;
opts.connect_list=[ones(nrepts,1),1+[1:nrepts]'];
opts.label_list='typenames';
%
opts=filldefault(opts,'color_norays','k');
opts=filldefault(opts,'LineWidth',1);
opts=filldefault(opts,'marker_origin','p');
opts=filldefault(opts,'color_origin','k');
opts=filldefault(opts,'axis_label_prefix','d');
%
sa=struct;
sa.typenames=labels;
%
%composite plot of mean response and lines to single-trials 
color_connect=cell(1,nrepts+1);
color_connect{1}=opts.color_norays;
for irept=1:nrepts
    color_connect{irept+1}=colors{irept};
end
opts.color_connect_sets_norays=color_connect;
opts.tag_text='includeTag';
opts_used=psg_plotcoords(cat(3,coords_mean(:,:),coords_trial(:,:,:)),dimlist,sa,struct(),opts);
%
%replot single trials to color the points
opts_trial=opts;
opts_trial=rmfield(opts_trial,'connect_list');
opts_trial=rmfield(opts_trial,'label_list');
opts_trial.axis_handle=opts_used.axis_handle;
opts_trial.tag_text='excludeTag'; %these don't go into the legends
for irept=1:nrepts
    opts_trial.color_norays=opts.color_connect_sets_norays{1+irept};
    opts_trial_used{irept}=psg_plotcoords(coords_trial(:,:,irept),dimlist,sa,struct(),opts_trial);
end
%
hp=plot3(0,0,0,opts.marker_origin); %plot origin
set(hp,'Color',opts.color_origin);
return
end

