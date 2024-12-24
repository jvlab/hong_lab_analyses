function hlid_rastim_trial_vis_legend%
% function hlid_rastim_trial_vis_legend
%
%utility legend creation for hlid_rastim_trial_vis
%
% note that this must be called after all plotting has been done, since it
% makes use of custom tags on handles (set up in hlid_rasti_trial_vis_plot)
% 
%   See also: HLID_RASTIM_TRIAL_VIS, HLID_RASTIM_TRIAL_VIS_PLOT.
%
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
for k=1:length(s1p1)
    %change ' set  1 connect   * point   1' to 'trial *'
    strings{k}=strrep(strrep(strrep(strrep(strrep(tags{s1p1(k)},'includeTag',' '),'point   1',''),'set  1',''),' ',''),'connect','trial ');
    handles(1,k)=hc(s1p1(k));
end
legend(handles,strings);
return
end
