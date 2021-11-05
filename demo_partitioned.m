% This scripts runs a benchmark to evaluate the performance of
% Multisync vs. edge averaging vs. full sync for partitioned 
% synchronization using the Roman Forum data set from 1DSfM
clc
clear
close all
warning off

addpath('./utils')
addpath('./multigraph')
addpath('./distributed_synch_v2')
addpath('./spectral_clustering')

% Fix RNG
rng('default');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.debug = true;
opts.max_shared = 50;
opts.data_name = 'Roman_Forum';
opts.test_runs = 1;
fprintf("Testing on dataset: %s\n\n", opts.data_name);

% general
opts.patchwork = true;

% SO(3) synchronization
opts.irls_global = false;
opts.irls = true;
opts.max_irls_iter = 30;
opts.tol=1e-5;
opts.weight_function='cauchy';

% clustering
opts.max_size_coefficient = 3.7;
opts.initial_datum = true;
opts.min_legs = 2;

% CEMP default params
opts.CEMP_parameters.max_iter = 100;
opts.CEMP_parameters.reweighting = 2.^((1:6)-1);
opts.CEMP_parameters.nsample = 50;

% MPLS default params
opts.MPLS_parameters.stop_threshold = 1e-3;
opts.MPLS_parameters.max_iter = 100;
opts.MPLS_parameters.reweighting = opts.CEMP_parameters.reweighting(end);
opts.MPLS_parameters.thresholding = [0.95,0.9,0.85,0.8];
opts.MPLS_parameters.cycle_info_ratio = 1./((1:opts.MPLS_parameters.max_iter)+1);

% Load dataset
[Z,A,R_gt,cameras_gt,i_cc,i_gt] = load_dataset(opts);


%% Partitioned synchronization
S_mean = [];
S_median = [];
M_mean = [];
M_median = [];

for i = 1:opts.test_runs
    results = run_snavely(Z,A,R_gt,cameras_gt,i_cc,i_gt,opts);
    S_mean = [S_mean; mean(results.err_standard)];
    S_median = [S_median; median(results.err_standard)];
    M_mean = [M_mean; mean(results.err_multi)];
    M_median = [M_median; median(results.err_multi)];
end

fprintf("\nTest results\n");
fprintf("Edge averaging\nMean: %.4f\nMedian: %.4f\n\n", mean(S_mean), mean(S_median));
fprintf("Multisync\nMean: %.4f\nMedian: %.4f\n\n", mean(M_mean), mean(M_median));