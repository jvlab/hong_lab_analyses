%hlid_kc_mdlcomp: Compare KC control and TNT, and Nest vs NoNest, for the
%16 stimuli in common with ORN data
% 
%   See also:  PSG_GET_COORDSETS, HLID_SETUP, HLID_ORN_KC_CHECK HLID_ORN_KC_MDLSUM.
%
%KC
hlid_setup;
opts_read.if_auto=1;
opts_read.if_log=0;
dmax=7;
%
stims_in_common=display_orders.kcmerge(1:end-1); %all but pfo
%
prefix_kc='./data/kc_tnt/hlid_odor17_coords_TNT*_consensus-pc_noscale-';
suffix_kc_nest='ovlp16sel-alpha.mat';
suffix_kc_nonest='ovlp16NoNest-alpha.mat';
%
disp('nested')
disp('control dataset')
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},{{cat(2,strrep(prefix_kc,'*','3c'),suffix_kc_nest)}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
kc_nest_3c=ds{1};
%
disp('label dataset')
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},{{cat(2,strrep(prefix_kc,'*','-label'),suffix_kc_nest)}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
kc_nest_label=ds{1};
%
disp('not nested')
disp('control dataset')
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},{{cat(2,strrep(prefix_kc,'*','3c'),suffix_kc_nonest)}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
kc_nonest_3c=ds{1};
%
disp('label dataset')
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},{{cat(2,strrep(prefix_kc,'*','-label'),suffix_kc_nonest)}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
kc_nonest_label=ds{1};
%
%Procrustes comparisons
disp('            3c vs label           Nest vs NoNest');
disp('      dim      Nest    NoNest      3c     label');
for k=1:dmax
    disp([k procrustes(kc_nest_3c{k},kc_nest_label{k}) procrustes(kc_nonest_3c{k},kc_nonest_label{k})...
        procrustes(kc_nest_3c{k},kc_nonest_3c{k}) procrustes(kc_nest_label{k},kc_nonest_label{k})]);
end


