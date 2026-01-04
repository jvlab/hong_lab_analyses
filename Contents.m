% Demos, pilots, etc for analysis of Hong Lab imaging data by J. Victor
%
%  File reading and preprocessing
%
% hlid_csv2coords_demo: convert model data in csv file to coordinates
% hlid_da_stimselect: select a stimulus set from a raw data file with several panels
% hlid_localopts: set up local options for file  names, etc
% hlid_nan_stims: script to set coords for some stimuli to NaN (mostly for debugging)
% hlid_orn_merge: merge ORN datasets keeping track of responses in each glomeruli
% hlid_rastim_trial_read: script for reading single-trial raw data
% hlid_rastim2coords_demo: demo conversion of response amplitudes to coordinates
% hlid_rastim2coords_pool: pool and convert response amplitudes to coordinates
% hlid_read_coorddata_demo: demo reading of a coordinate file
% hlid_remove_stims: script to remove stimuli from a coordinate file (mostly for debugging)
% hlid_setup: set up options for Hong Lab imaging data, to run psg_visualize_demo, psg_consensus_demo, _psg_coord_pipe_proc
% hlid_setup_func: function called by hlid_setup
%
%  Analysis
%
% hlid_connect_eigstats: statistics of eigenvalues of connectivity matrices
% hlid_connect_eigstats2: comparison of eigenvectors of connectivity matrices
% hlid_csv2connectivity: read connectivity data
% hlid_csv2connectivity_demo: examine connectivity data
% hlid_coords_plot: utility to plot coords and SVD after conversion from data files
% hlid_coords_svd: utility to create coordinates from data files via SVD, and to create associated metadata
% hlid_geom_transform_stats_label: utility labeling for hlid_geom_traansform_stats_summ
% hlid_geom_transform_stats: analysis of transformations between datasets, with shuffle statistics
% hlid_geom_transform_stats_3dplot: plot of representational spaces for hlid_geom_transform_stats
% hlid_geom_transform_stats_summ: summary plot for hli_geom_transform_stats
% hlid_kc_mdlcomp: compare kc space constructed via nested or non-nested pipeline
% hlid_majaxes: examine axes identified in a transformation
% hlid_majaxes_compare: compare axes identified in a transformation with axes from connectivity
% hlid_majaxes_eval:  compare axes identified in a transformation with principal components
% hlid_majaxes_eval2: as in hlid_majazxes_eval, but add analysis of fitted transformation, other params
% hlid_mds_coords_geomodels: transformations between spaces constructed via mds and nonstandard embeddigs, e.g., Pearson, cosine
% hlid_mds_transform_jackstats: statistics related to hlid_mds_coords_geomodels, jackknife by stimuli
% hlid_mds_transform_jackstats_summ: summary and plots of hlid_mds_transform_jackstats: magnif factors
% hlid_mds_transform_jackstats_summ2: summary and plots of hlid_mds_transform_jackstats: direction of axes of transformation
% hlid_mds_transform_stats: statistics related to hlid_mds_coords_geomodels: shuffles and bootstraps
% hlid_orn_kc_check: check on variances, etc, from coordinates for orn and kc for overlapping stimuli
% hlid_orn_kc_mdlsum: basic goodness of fit via Procrustes for orn->kc transformation
% hlid_participation_ratio: compute participation ratio from single-trial datasets
% hlid_pool_pcacontrib: show contributions of each component dataset to pcs of a pooled dataset
% hlid_rastim_mds_coords_demo: examine nonstandard embeddings, e.g., cosine, Pearson, via transformations of z-scores, prior to computing coordinates
% hlid_rastim_mds_coords_subsamp: calculations of hlid_rastim_coords_demo on subsamples of files
% hlid_rastim_mds_coords_summ: summarize output statistics of hlid_rastim_mds_coords_demo
% hlid_rastim_mds_coords_summ2: unexplained variance in hlid_rastim_mds_coords_demo: compare orn and kc
% hlid_rastim_mds_coords_summ2s: plots results of hlid_rastim_mds_coords_subsamp
% hlid_rastim_trial_decode: single-trial decoding of raw data
% hlid_rastim_trial_decode_confmtx: confusion matrix plotting for hlid_rastim_trial_decode
% hlid_rastim_trial_decode_3dplot: resonse space plotting for hlid_rastim_trial_decode
% hlid_rastim_trial_pca: single-trial analysis of raw data
% hlid_rastim_trial_pca2: further plots for single-trial analysis of raw data
% hlid_rastim_trial_plot: utility plotting for hlid_rastim_trial_pca
% hlid_rastim_trial_vis: single-trial visualization of raw data
% hlid_rastim_trial_vis2: further analysis (stats of covariances) for hlid_rastim_trial_vis
% hlid_rastim_trial_vis_legend: plotting utility plotting for hlid_rastim_trial_vis
% hlid_rastim_trial_vis_plot: plotting utility plotting for hlid_rastim_trial_vis
% hlid_tnt_label_model_3c_compare: compare control kc dataset with tnt and tnt transformed by models
%
%   Copyright (c) 2024, 2025, 2026 by J. Victor
