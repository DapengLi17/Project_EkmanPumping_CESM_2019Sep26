clear all;close all;clc;

addpath(genpath('eddy_detect'));
% addpath(genpath('/scratch/user/sanjiv/Matlab-code/matlabtools/'));
% addpath(genpath('/scratch/user/sanjiv/Matlab-code/ROMS_Jaison/'));
% addpath(genpath('/scratch/user/sanjiv/Matlab-code/matlabtools/'));

gridfile = '../data_after_manipulation/POP2ROMS_mask_2020Aug05.nc';
% ncgr = netcdf(gridfile,'read')
% grid.lat      = ncgr{'lat_rho'}(:) ;
reffile = '../data_after_manipulation/POP2ROMS_Monthly_TAUXTAUYSSH_90-01.nc';
omod       = 'basic';

grid = func_init_grid (gridfile, reffile, omod)