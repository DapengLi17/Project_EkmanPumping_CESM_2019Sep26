%% ===readme===

% descrip: matlab scripts plot Nonlinear Ekman vertical heat flux Q
% based on Stern (1965) formulation  
% 2) Stern nonlinear Ekman vertical pumping velocity
% 1) linear Ekman vertical pumping velocity
% 3) CESM model output

% update history:
% v1.0 DL 2019Dec05
% v1.1 DL 2019Dec09
% v1.2 DL 2020Feb05

% extra notes:
% =============


%% === set up environments ===
clear all;close all;clc;

  date_str='2020Feb11'; 

% --- summer seasons ---
header = 'This file was generated by DL via code ComputeSaveQekmanEddyGlobalBorealSummerAv_2020Feb11.m';

MonSea = [4:9];
infileWceAv = ['../data_after_manipulation/CESM_regrid_Apr-SepAv_WVEL_Global.nc'];
infileWekAv = ['../data_after_manipulation/TauKesaiRoWvelTempGlobalBorealSummerAv_2020Feb08.mat'];
infileTwAv  = ['../data_after_manipulation/CESM_regrid_Apr-SepAv_TEMP_Global.nc'];

outfileQ = ['../data_after_manipulation/QekmanEddyGlobalBorealSummerAv_',date_str,'.mat'];
% ----------------------

addpath(genpath('Func4EkmanProject/'))

Constants4CESM_GlobalDpthProf

wt_nl_prim_csum = 0;
wt_li_prim_csum = 0;
wt_ce_prim_csum = 0;
Ndays = 0; % number of loops
rho_w = 1020;
Cp    = 3.99e3; % [J/kg/k] the specific heat at constant pressure (Thorpe 2007)
%=============================


%% === data loading and analysis ===
EkAv = load(infileWekAv);
lon_1d  = EkAv.lon_1d;
lat_1d  = EkAv.lat_1d;
w_nl_av = EkAv.w_nl_TimeAv;
w_li_av = EkAv.w_li_TimeAv;
clear EkAv

w_ce_av = ncread(infileWceAv,'WVEL',start_4dvar,count_4dvar,stride_4dvar);
z_wvel_ce = double(ncread(infileWvelNc,'z_w_top',start_z,count_z,stride_z)).*0.01;

Tw_av   = ncread(infileTwAv,'TEMP',start_4dvar,count_4dvar,stride_4dvar);
z_Tw = double(ncread(infileTwAv,'z_t',start_z,count_z,stride_z)).*0.01;

% read lon, lat, and z (non-change dat during loop)

for iyr = 87:90
  for imon = MonSea
      
%% === load data ===
disp(['working on ',num2str(iyr),'-',num2str(imon,'%02d')])

infileAvMat  = ['../data_after_manipulation/WvelTempGlobal_CESM_', ...
    num2str(iyr),'-',num2str(imon,'%02d'),'.mat'];
Ek = load(infileAvMat);
w_nl_3d = Ek.w_nl_3d;
w_li_3d = Ek.w_li_3d;
clear infileMat Ek


infileWvelNc  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/WVEL/' ...
    'CESM_regrid_MonthlyOutput_WVEL_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
infileTempNc  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/' ...
    'CESM_regrid_MonthlyOutput_TEMP_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)

wvel = ncread(infileWvelNc,'WVEL', start_4dvar,count_4dvar,stride_4dvar).*0.01;
Tw   = ncread(infileTempNc,'TEMP', start_4dvar,count_4dvar,stride_4dvar);
% time = ncread(infile_TTS, 'time');
ntime = size(w_nl_3d,3); % time dimension
clear infileWvelNc infileTempNc
% ==================

% loop through time for .mat files
for iday = 1 : ntime 
disp(['analyzing day ', num2str(iday,'%02d')])

% --- compute w', T', w'T' --- 
w_nl_prim = w_nl_3d(:,:,iday)-w_nl_av;
w_li_prim = w_li_3d(:,:,iday)-w_li_av;

w_ce_prim = wvel(:,:,:,iday) - w_ce_av;
Tw_prim = Tw(:,:,:,iday) - Tw_av;

wt_nl_prim = w_nl_prim.*Tw_prim;
wt_li_prim = w_li_prim.*Tw_prim;
wt_ce_prim = w_ce_prim.*Tw_prim;
% -----------------------------

% --- compute cumulative sum (csum) --- 
wt_nl_prim_csum = wt_nl_prim_csum + wt_nl_prim;
wt_li_prim_csum = wt_li_prim_csum + wt_li_prim;
wt_ce_prim_csum = wt_ce_prim_csum + wt_ce_prim;
% ---------------

% --- compute number of loops ---
Ndays = Ndays + 1;
% -------------------------------



%% === data analysis ===
for iday = 1 : ntime % loop through time
disp(['analyzing day ', num2str(iday,'%02d')])

% --- compute w', T', w'T' --- 
w_nl_prim = w_nl_4d(:,:,iday)-w_nl_av;
w_li_prim = w_li_4d(:,:,iday)-w_li_av;
w_ce_prim = w_ce_4d(:,:,iday)-w_ce_av;
Tw_prim   = Tw_3d(:,:,iday)-Tw_av;

wt_nl_prim = w_nl_prim.*Tw_prim;
wt_li_prim = w_li_prim.*Tw_prim;
wt_ce_prim = w_ce_prim.*Tw_prim;
% -----------------------------

% --- compute cumulative sum (csum) --- 
wt_nl_prim_csum = wt_nl_prim_csum + wt_nl_prim;
wt_li_prim_csum = wt_li_prim_csum + wt_li_prim;
wt_ce_prim_csum = wt_ce_prim_csum + wt_ce_prim;
% ---------------

% --- compute number of loops ---
Ndays = Ndays + 1;
% -------------------------------

% --- clear vars to save memory --- 
clear w_nl_prim w_li_prim w_ce_prim Tw_prim 
clear wt_nl_prim wt_li_prim wt_ce_prim 
% ---------------------------------
end

% --- clear vars to save memory --- 
  clear w_nl_3d w_li_3d w_ce_3d Tw_3d ntime 
% ---------------------------------

  end
end

% --- compute mean vertical heat flux ---
Q_nl = rho_w.*Cp.*wt_nl_prim_csum ./ Ndays;
Q_li = rho_w.*Cp.*wt_li_prim_csum ./ Ndays;
Q_ce = rho_w.*Cp.*wt_ce_prim_csum ./ Ndays;
% --------------------------
% ==================
  
  
%% === output data ===
save(outfileQ,'header','lon_1d','lat_1d','Ndays', ...
  'z_wvel_ce','z_Tw','Q_nl','Q_li','Q_ce') 
% ====================
