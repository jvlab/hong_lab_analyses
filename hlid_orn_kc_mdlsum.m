%hlid_orn_kc_mdlsum: Summarize modeling of orn->kc transformation
%
% Note that model is affine with offset, fitted to centered data.
% On the kc side, data are centered when finding a consensus across flies and
% these data files are used to fit the model, then pfo is removed so the data in the
% *TNT-label_consensus-pc_noscale-ovlp16sel-alpha files are not centered.
% 
% On the orn side, data are centered internally in psg_geomodels_run based on the 16 stimuli in common
% and then recentered when the model is applied.
%
% For goodness of fit, centering does not matter as measured by Procrustes.
%
%   See also:  PSG_MAJAXES_EVAL, PSG_GEOMODELS_RUN, PSG_MAJAXES_REORDER PSG_GET_COORDSETS, HLID_SETUP, HLID_ORN_KC_CHECK.
%
%
hlid_setup;
opts_read.if_auto=1;
opts_read.if_log=0;
dmax=7;
prefix_orn='./data/orn_merged/hlid_odor039_coords_merged_union-ovlp16sel-alpha';
prefix_kc='./data/kc_tnt/hlid_odor17_coords_TNT*_consensus-pc_noscale-';
suffix_kc='ovlp16sel-alpha.mat';
%
stims_in_common=display_orders.kcmerge(1:end-1); %all but pfo
%
disp(' ');
disp('****************');
disp('ORN datasets used to assess modeling');
disp('merged');
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{prefix_orn}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_mdlsum_util(ds{1},dmax);
orn_merged=ds{1};
%
disp('transformed by orn->TNT3c model');
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{cat(2,prefix_orn,'_mdl3c')}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_mdlsum_util(ds{1},dmax);
orn_model_3c=ds{1};
%
disp('transformed by orn->TNTlabel model');
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{cat(2,prefix_orn,'_mdl-label')}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_mdlsum_util(ds{1},dmax);
orn_model_label=ds{1};
%
disp('****************');
disp('KC datasets used to assess modeling');
disp('control');
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{cat(2,strrep(prefix_kc,'*','3c'),suffix_kc)}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_mdlsum_util(ds{1},dmax);
kc_3c=ds{1};
%
disp('label');
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{cat(2,strrep(prefix_kc,'*','-label'),suffix_kc)}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_mdlsum_util(ds{1},dmax);
kc_label=ds{1};

disp('********************')
disp('procrustes evaluation of fits');
disp(' orn to control kc')
disp(' dim         merged    merged+xform');
for k=1:dmax
    disp([k procrustes(kc_3c{k},orn_merged{k}) procrustes(kc_3c{k},orn_model_3c{k})]);
end
disp(' orn to TNT-label kc')
disp(' dim         merged    merged+xform');
for k=1:dmax
    disp([k procrustes(kc_label{k},orn_merged{k}) procrustes(kc_label{k},orn_model_label{k})]);
end

function hlid_orn_kc_mdlsum_util(ds,dtest)

%check nesting
c=ds{dtest};
cfull=ds{end};
nest_check=max(max(abs(c-cfull(:,1:dtest))));
disp(sprintf(' nest check: %15.10f',nest_check));
means=mean(c,1);
disp(sprintf(' mean check: %15.10f',max(abs(means))));
% disp('means');
% disp(sprintf(' %15.10f',means));
% vars_zero=mean(c.^2,1);
% disp('variances around zero and square roots');
% disp(sprintf(' %15.10f',vars_zero));
% disp(sprintf(' %15.10f',sqrt(vars_zero)));
orth_check_zero=max(max(abs(c'*c-diag(diag(c'*c)))));
disp(sprintf(' orthog check around zero: %15.10f',orth_check_zero));
% vars_ms=mean(c_ms.^2,1);
% disp('mean-subtracted variances and square roots');
% disp(sprintf(' %15.10f',vars_ms));
% disp(sprintf(' %15.10f',sqrt(vars_ms)));
c_ms=c-repmat(means,size(c,1),1);
orth_check_ms=max(max(abs(c_ms'*c_ms-diag(diag(c_ms'*c_ms)))));
disp(sprintf(' mean-subtracted orthog check: %15.10f',orth_check_ms));
return
end
