%hlid_majaxes: analyze the transformations of an affine geometric model to determine major axes
% for Hong Lab data
%
% needs two data files (one as reference, one to be adjusted), and the results of psg_geomodels_run that
% determined geometric models between adjusted and reference dataset
%
%   See also: PSG_GEOMODELS_RUN, PSG_GEOMODELS_DEFINE, HLID_SETUP, PSG_MAJAXES.
%
hlid_setup;
%
if ~exist('opts_read') opts_read=struct;end
opts_read=filldefault(opts_read,'if_log',1);
opts_read=filldefault(opts_read,'if_justsetup',0);
opts_read=filldefault(opts_read,'data_fullname_def',hlid_opts.coord_data_fullname_def);
opts_read=filldefault(opts_read,'setup_fullname_def',[]);
opts_read=filldefault(opts_read,'domain_list',hlid_opts.domain_list_def);
opts_read=filldefault(opts_read,'setup_suffix',hlid_opts.setup_setup_suffix);
opts_read=filldefault(opts_read,'type_class_aux',hlid_opts.type_class_def);
%
if ~exist('opts_majaxes') opts_majaxes=struct;end
opts_majaxes=filldefault(opts_majaxes,'plot_pairs',[3 3;4 3;5 3]);
opts_majaxes=filldefault(opts_majaxes,'plot_order',display_orders.kcclust);
%
if ~exist('fn_ref') fn_ref='./data/kc_soma_nls/hlid_odor17_coords_kc_soma_nls_consensus_scale.mat'; end
if ~exist('fn_adj') fn_adj='./data/orn_terminals/hlid_odor17_coords_orn_terminals_consensus_scale.mat'; end
if ~exist('fn_geo') fn_geo='psg_geomodels_run_orn7setsConsensusScale_kcConsensusScale_6d.mat'; end
% alternate: psg_geomodels_run_orn7setsPooled_kcPooled.mat % this has piecewise affine also but is pooled, not consensus
%
fn_ref=getinp('file name for dataset that is the reference','s',[],fn_ref);
[d_ref,sa_ref]=psg_read_coorddata(fn_ref,[],opts_read);
%
fn_adj=getinp('file name for dataset that is adjusted','s',[],fn_adj);
ds_adj=load(fn_adj);
[d_adj,sa_adj]=psg_read_coorddata(fn_adj,[],opts_read);
%
fn_geo=getinp('file name containing results of psg_geomodels_run','s',[],fn_geo);
results_geo=getfield(load(fn_geo),'results');
%
[results_axes,opts_majaxes_used]=psg_majaxes(d_ref,sa_ref,d_adj,sa_adj,results_geo,opts_majaxes);
