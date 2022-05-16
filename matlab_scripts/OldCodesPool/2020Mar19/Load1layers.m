clear all;close all;clc;
infile1 = '../raw_data/CESMncFilesGlobalFromAgartha_ZhangQiuying/CESM_regrid_MonthlyOutput_WVEL_87-03.nc';
% ncdisp(infile1)

infile2 = '../raw_data/CESMncFilesGlobalFromAgartha_ZhangQiuying/CESM_regrid_MonthlyOutput_WVEL_87-03_z_w_top_equap_5.nc';
% ncdisp(infile2)

Constants4CESM_Global

wvel1 = squeeze(ncread(infile1,'WVEL',start_4dvar,count_4dvar,stride_4dvar)).*0.01;
z1     = ncread(infile1,'z_w_top',start_z,count_z,stride_z)
wvel2 = squeeze(ncread(infile2,'WVEL')).*0.01;
z2 = squeeze(ncread(infile2,'z_w_top'))

dif_wvel = wvel2 - wvel1;
a = dif_wvel(:);
nansum(abs(a))