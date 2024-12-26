function hlid_rastim_trial_vis_legend(maxtrial)
% function hlid_rastim_trial_vis_legend(maxtrial)
%
%utility legend creation for hlid_rastim_trial_vis
%
% maxtrial, if supplied, is largest trial number to show in legend
%
% note that this must be called after all plotting has been done, since it
% makes use of custom tags on handles (set up in hlid_rasti_trial_vis_plot)
% 
%   See also: HLID_RASTIM_TRIAL_VIS, HLID_RASTIM_TRIAL_VIS_PLOT.
%
if (nargin<1)
    maxtrial=Inf;
end
hc=get(gca,'Children');
tags=cell(length(hc),1);
for ich=1:length(hc)
    tags{ich}=get(hc(ich),'Tag');
end
hc_incl=find(contains(tags,'includeTag'));
hc_p1=find(contains(tags,'point   1'));
hc_s1=find(contains(tags,'set  1'));
hc_exclude=find(contains(tags,'excludeTag'));
s1p1=setdiff(intersect(hc_incl,intersect(hc_p1,hc_s1)),hc_exclude);
handles=[];
strings=cell(0);
trialnos=zeros(1,length(s1p1));
for k=1:length(s1p1)
    %change ' set  1 connect   * point   1' to 'trial *'
    strings{k,1}=strrep(strrep(strrep(strrep(strrep(tags{s1p1(k)},'includeTag',' '),'point   1',''),'set  1',''),' ',''),'connect','trial ');
    trialnos(k)=str2num(strrep(strings{k,1},'trial',''));
    handles(k,1)=hc(s1p1(k));
end
legend(flipud(handles(trialnos<=maxtrial)),flipud(strings(trialnos<=maxtrial)));
return
end
