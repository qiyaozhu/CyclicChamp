%% This script is for testing the Backbone_SA.m
% Test 1: sample backbones without disulfide bond
clc;
clear all;

n = 12; % number of backbone residues;
N = 50; % number of simulated annealing trajectories to try
best_k = 1; % for each trajectory, write out the best_k backbones found as pdbs
outdir = "Backbones_nodisulfide/"; % the output directory

is_disulfide = false;
no_cys = [];
min_sep = -1;
max_sep = -1;

Backbone_SA(n, N, best_k, outdir, is_disulfide, no_cys, min_sep, max_sep);


%% Test 2: sample backbones with disulfide bond
clc;
clear all;

n = 12; % number of backbone residues;
N = 50; % number of simulated annealing trajectories to try
best_k = 1; % for each trajectory, write out the best_k backbones found as pdbs
outdir = "Backbones_disulfide/"; % the output directory

is_disulfide = true;
no_cys = [1, 5]; % positions to be excluded from forming disulfide bond
min_sep = floor(n/4)+1; % min separation between the two cys residues, bond (i, i+sep)
max_sep = floor(n/2);

Backbone_SA(n, N, best_k, outdir, is_disulfide, no_cys, min_sep, max_sep);

