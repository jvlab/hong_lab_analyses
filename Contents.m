% Demos, pilots, etc for analysis of Hong Lab imaging data by J. Victor
%
% hlid_connect_eigstats: statistics of eigenvalues of connectivity matrices
% hlid_connect_eigstats2: comparison of eigenvectors of connectivity matrices
% hlid_csv2connectivity: read connectivity data
% hlid_csv2connectivity_demo: examine connectivity data
% hlid_csv2coords_demo: convert model data in csv file to coordinates
% hlid_coords_plot: utility to plot coords and SVD after conversion from data files
% hlid_coords_svd: utility to create coordinates from data files via SVD, and to create associated metadata
% hlid_localopts: set up local options for file  names, etc
% hlid_majaxes: examine axes identified in a transformation
% hlid_majaxes_compare: compare axes identified in a transformation with axes from connectivity
% hlid_nan_stims: script to set coords for some stimuli to NaN (mostly for debugging)
% hlid_orn_merge: merge ORN datasets keeping track of responses in each glomeruli
% hlid_pool_pcacontrib: show contributions of each component dataset to pcs of a pooled dataset
% hlid_rastim_trial_pca: single-trial analysis of raw data
% hlid_rastim_trial_pca2: further plots for single-trial analysis of raw data
% hlid_rastim_trial_plot: utility plotting for hlid_rastim_trial_pca
% hlid_rastim_trial_read: script for reading single-trial raw data
% hlid_rastim_trial_vis: single-trial visualization of raw data
% hlid_rastim_trial_vis2: further analysis (stats of covariances) for hlid_rastim_trial_vis
% hlid_rastim_trial_vis_legend: plotting utility plotting for hlid_rastim_trial_vis
% hlid_rastim_trial_vis_plot: plotting utility plotting for hlid_rastim_trial_vis
% hlid_rastim2coords_demo: demo conversion of response amplitudes to coordinates
% hlid_rastim2coords_pool: pool and convert response amplitudes to coordinates
% hlid_read_coorddata_demo: demo reading of a coordinate file
% hlid_remove_stims: script to remove stimuli from a coordinate file (mostly for debugging)
% hlid_setup: set up options for Hong Lab imaging data, to run psg_visualize_demo, psg_consensus_demo, _psg_coord_pipe_proc
%
%   Copyright (c) 2024, 2025 by J. Victor
