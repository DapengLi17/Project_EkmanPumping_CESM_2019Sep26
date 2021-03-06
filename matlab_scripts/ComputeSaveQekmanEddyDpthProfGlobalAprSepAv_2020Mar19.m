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

  date_str='2020Mar19'; 

% --- summer seasons ---
header = 'This file was generated by DL via code ComputeSaveQekmanEddyDpthProfGlobalAprSepAv_2020Mar19.m';

MonSea = [4:9]; 
infileWekAv = ['../data_after_manipulation/TauKesaiRoWvelTempGlobalAprSepAv_2020Mar17.mat'];
infileWceAv = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/WVEL/CESM_regrid_Apr-SepAv_WVEL_Global.nc'];
infileTwAv  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/CESM_regrid_Apr-SepAv_TEMP_Global.nc'];

outfileQ = ['../data_after_manipulation/QekmanEddyDpthProfGlobalAprSepAv_',date_str,'.mat'];
% ----------------------

addpath(genpath('Func4EkmanProject/'))

Constants4CESM_GlobalDpthProf

wt_te_prim_csum = 0;
wt_ne_prim_csum = 0;
wt_le_prim_csum = 0;
wt_ce_prim_csum = 0;
Ndays = 0; % number of loops
rho_w = 1020;
Cp    = 3.99e3; % [J/kg/k] the specific heat at constant pressure (Thorpe 2007)
%=============================


%% === data loading and analysis ===
EkAv = load(infileWekAv);
lon_1d  = EkAv.lon_1d;
lat_1d  = EkAv.lat_1d;
w_te_av = EkAv.w_te_TimeAv;
w_ne_av = EkAv.w_ne_TimeAv;
w_le_av = EkAv.w_le_TimeAv;
clear EkAv

w_ce_av = ncread(infileWceAv,'WVEL',start_4dvar,count_4dvar,stride_4dvar).*0.01;
z_wvel_ce = double(ncread(infileWceAv,'z_w_top',start_z,count_z,stride_z)).*0.01;

Tw_av   = ncread(infileTwAv,'TEMP',start_4dvar,count_4dvar,stride_4dvar);
z_Tw = double(ncread(infileTwAv,'z_t',start_z,count_z,stride_z)).*0.01;

% read lon, lat, and z (non-change dat during loop)

for iyr = 87:90
  for imon = MonSea
      
%% === load data ===
disp(['working on ',num2str(iyr),'-',num2str(imon,'%02d')])

infileMat  = ['../data_after_manipulation/WvelTempGlobal_CESM_', ...
    num2str(iyr),'-',num2str(imon,'%02d'),'.mat'];
Ek = load(infileMat);
w_te_3d = Ek.w_te_3d;
w_ne_3d = Ek.w_ne_3d;
w_le_3d = Ek.w_le_3d;
clear infileMat Ek

infileWvelNc  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/WVEL/' ...
    'CESM_regrid_MonthlyOutput_WVEL_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
infileTempNc  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/' ...
    'CESM_regrid_MonthlyOutput_TEMP_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)

w_ce_4d = ncread(infileWvelNc,'WVEL', start_4dvar,count_4dvar,stride_4dvar).*0.01;
Tw_4d   = ncread(infileTempNc,'TEMP', start_4dvar,count_4dvar,stride_4dvar);
% time = ncread(infile_TTS, 'time');
ntime = size(w_ne_3d,3); % time dimension
clear infileWvelNc infileTempNc
% ==================

% loop through time for .mat files
for iday = 1 : ntime 
disp(['analyzing day ', num2str(iday,'%02d')])

% --- compute w', T', w'T' --- 
w_te_prim = [w_te_3d(:,:,iday)-w_te_av]'; % size: 3600x1401
w_ne_prim = [w_ne_3d(:,:,iday)-w_ne_av]';
w_le_prim = [w_le_3d(:,:,iday)-w_le_av]';

w_te_prim_rep = repmat(w_te_prim,[1,1,count_z]); % size: 3600x1401x9
w_ne_prim_rep = repmat(w_ne_prim,[1,1,count_z]);
w_le_prim_rep = repmat(w_le_prim,[1,1,count_z]);

w_ce_prim = w_ce_4d(:,:,:,iday) - w_ce_av; % size: 3600x1401x9 
Tw_prim = Tw_4d(:,:,:,iday) - Tw_av;

wt_te_prim = w_te_prim_rep.*Tw_prim;
wt_ne_prim = w_ne_prim_rep.*Tw_prim;
wt_le_prim = w_le_prim_rep.*Tw_prim;
wt_ce_prim = w_ce_prim.*Tw_prim;
% -----------------------------

% --- compute cumulative sum (csum) --- 
wt_te_prim_csum = wt_te_prim_csum + wt_te_prim;
wt_ne_prim_csum = wt_ne_prim_csum + wt_ne_prim;
wt_le_prim_csum = wt_le_prim_csum + wt_le_prim;
wt_ce_prim_csum = wt_ce_prim_csum + wt_ce_prim;
% -------------------------------------

% --- compute number of loops ---
Ndays = Ndays + 1;
% -------------------------------

% --- clear vars to save memory --- 
clear w_te_prim w_ne_prim w_le_prim 
clear w_te_prim_rep w_ne_prim_rep w_le_prim_rep w_ce_prim Tw_prim 
clear wt_te_prim wt_ne_prim wt_le_prim wt_ce_prim 
% ---------------------------------
end

% --- clear vars to save memory --- 
  clear w_te_3d w_ne_3d w_le_3d w_ce_4d Tw_4d ntime 
% ---------------------------------

  end
end

% --- compute mean vertical heat flux ---
Q_te = rho_w.*Cp.*wt_te_prim_csum ./ Ndays;
Q_ne = rho_w.*Cp.*wt_ne_prim_csum ./ Ndays;
Q_le = rho_w.*Cp.*wt_le_prim_csum ./ Ndays;
Q_ce = rho_w.*Cp.*wt_ce_prim_csum ./ Ndays;
% --------------------------
% ==================
  
  
%% === output data ===
save(outfileQ,'header','lon_1d','lat_1d','Ndays', ...
  'z_wvel_ce','z_Tw','Q_te','Q_ne','Q_le','Q_ce') 
% ====================
