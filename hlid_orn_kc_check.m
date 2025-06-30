%hlid_orn_kc_check: Do some basic checks on the merged ORN files and KC-TNT files,
%for comparison with outputs of hlid_majaxes_eval
%
% Note that hlid_majaxes_eval uses output of psg_geomodels_run, which centers the datasets around
% all stimuli, so the variances around zero computed by hlid_majaxes_eval will differ from
% the variances around zero computed here. However, for the 16 stimuli in
% common, the sqrt(variances) around the mean should match.
%
%   See also:  PSG_MAJAXES_EVAL, PSG_GEOMODELS_RUN, PSG_MAJAXES_REORDER PSG_GET_COORDSETS, HLID_SETUP.
%
if ~exist('dtest') dtest=4; end
hlid_setup;
opts_read.if_auto=1;
opts_read.if_log=0;
%
stims_in_common=display_orders.kcmerge(1:end-1); %all but pfo
%
disp(' ');
disp('****************');
disp('ORN datasets merged from megamat_wd and validation2_wd, and then reduced to 16 stimuli');
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{'./data/orn_merged/hlid_odor039_coords_merged_union.mat'}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_check_util(ds{1},dtest);
%
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{'./data/orn_merged/hlid_odor039_coords_merged_union-ovlp16.mat'}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_check_util(ds{1},dtest);
%
disp('restrict to the 16 stimuli common to orns and kc-tnt')
ds_restrict=cell(1,length(ds{1}));
for id=1:length(ds{1})
    [ds_restrict{id},xtick_labels]=psg_majaxes_reorder(ds{1}{id},stims_in_common,sas{1}.typenames);
end
hlid_orn_kc_check_util(ds_restrict,dtest);
%
disp(' ');
disp('****************');
disp('kc_tnt datasets: control: consensus and then reduced to 16 stimuli');
%
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{'./data/kc_tnt/hlid_odor17_coords_TNT3c_consensus-pc_noscale.mat'}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_check_util(ds{1},dtest);
%
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{'./data/kc_tnt/hlid_odor17_coords_TNT3c_consensus-pc_noscale-ovlp16.mat'}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_check_util(ds{1},dtest);
disp('restrict to the 16 stimuli common to orns and kc-tnt')
ds_restrict=cell(1,length(ds{1}));
for id=1:length(ds{1})
    [ds_restrict{id},xtick_labels]=psg_majaxes_reorder(ds{1}{id},stims_in_common,sas{1}.typenames);
end
hlid_orn_kc_check_util(ds_restrict,dtest);
%
disp(' ');
disp('****************');
disp('kc_tnt datasets: label: consensus and then reduced to 16 stimuli');
%
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{'./data/kc_tnt/hlid_odor17_coords_TNT-label_consensus-pc_noscale.mat'}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_check_util(ds{1},dtest);
%
[sets,ds,sas,rayss,opts_read_used]=psg_get_coordsets(setfields(opts_read,{'data_fullnames'},...
    {{'./data/kc_tnt/hlid_odor17_coords_TNT-label_consensus-pc_noscale-ovlp16.mat'}}),[],[],1);
disp(opts_read_used{1}.data_fullnames);
hlid_orn_kc_check_util(ds{1},dtest);
disp('restrict to the 16 stimuli common to orns and kc-tnt')
ds_restrict=cell(1,length(ds{1}));
for id=1:length(ds{1})
    [ds_restrict{id},xtick_labels]=psg_majaxes_reorder(ds{1}{id},stims_in_common,sas{1}.typenames);
end
hlid_orn_kc_check_util(ds_restrict,dtest);

function hlid_orn_kc_check_util(ds,dtest)

%check nesting
c=ds{dtest};
cfull=ds{end};
nest_check=max(max(abs(c-cfull(:,1:dtest))));
disp(sprintf(' nest check: %15.10f',nest_check));
means=mean(c,1);
disp('means');
disp(sprintf(' %15.10f',means));
vars_zero=mean(c.^2,1);
disp('variances around zero and square roots');
disp(sprintf(' %15.10f',vars_zero));
disp(sprintf(' %15.10f',sqrt(vars_zero)));
orth_check_zero=max(max(abs(c'*c-diag(diag(c'*c)))));
disp(sprintf(' orthog check around zero: %15.10f',orth_check_zero));
c_ms=c-repmat(means,size(c,1),1);
vars_ms=mean(c_ms.^2,1);
disp('mean-subtracted variances and square roots');
disp(sprintf(' %15.10f',vars_ms));
disp(sprintf(' %15.10f',sqrt(vars_ms)));
orth_check_ms=max(max(abs(c_ms'*c_ms-diag(diag(c_ms'*c_ms)))));
disp(sprintf(' mean-subtracted orthog check: %15.10f',orth_check_ms));
return
end
