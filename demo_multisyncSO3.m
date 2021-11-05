% This scripts runs a benchmark to evaluate the performance of 
% Multisync vs. edge averaging in SO(3)
clc
close all
clear

% Add script paths
addpath('./utils/')
addpath('./multigraph/')
addpath('./bench/')

% Fix RNG
rng('default')

% Test params
N = 10;
m_avg = 10;
sigma = pi/12;
runs = 10;


%% Evaluate performance

% Run the benchmark
S = [];
M = [];
for i = 1:runs
    fprintf("Run %d of %d\n", i, runs);
    [s, m] = simple_multi_bench(N, m_avg, sigma);
    S = [S; s'];
    M = [M; m'];
end

% Print the results
fprintf("\nTest results\n");
fprintf("Edge averaging:\n%f mean\n%f median\n%f var\n\n", ...
    mean(S), median(S), var(S));
fprintf("Multisync:\n%f mean\n%f median\n%f var\n\n", ...
    mean(M), median(M), var(M));