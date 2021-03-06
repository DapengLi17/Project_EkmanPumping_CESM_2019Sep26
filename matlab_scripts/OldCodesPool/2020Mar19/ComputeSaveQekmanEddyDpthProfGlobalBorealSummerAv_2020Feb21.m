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
infileWavProf  = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/CESM_regrid_Apr-SepAv_WVEL_Global.nc';
infileTavProf = '../data_after_manipulation/CESM_regrid_Oct-MarAv_TEMP_Global.nc';

outfileQ = ['../data_after_manipulation/QekmanEddyGlobalBorealSummerAv_',date_str,'.mat'];
MonSea = [4:9]; 
% ----------------------

addpath(genpath('Func4EkmanProject/'))

% Constants4CESM_Global

wt_nl_prim_csum = 0;
wt_li_prim_csum = 0;
wt_ce_prim_csum = 0;
Ndays = 0; % number of loops
rho_w = 1020;
Cp    = 3.99e3; % [J/kg/k] the specific heat at constant pressure (Thorpe 2007)
%=============================


%% === data loading and analysis ===
% --- load av fields ---
Av = load(infileAv);
w_nl_av = Av.w_nl_TimeAv;
w_li_av = Av.w_li_TimeAv;
lon_1d  = Av.lon_1d;
lat_1d  = Av.lat_1d;

clear Av

z_wel_ce= ncread(infileWav,'z_');
z_Tw    = ncread(infileTav,'z_');
w_ce_av = ncread(infileWav,'WVEL');
Tw_av   = ncread(infileTavProf,'TEMP');
% ----------------------

for iyr = 87:90
  for imon = MonSea
      
%% === load data ===
infileMat  = ['../data_after_manipulation/WvelTempGlobal_CESM_', ...
    num2str(iyr),'-',num2str(imon,'%02d'),'.mat'];
disp(['loading ', infileMat])
dummy = load(infileMat);
w_nl_3d = dummy.w_nl_3d;
w_li_3d = dummy.w_li_3d;
clear infile dummy
ntime = size(Tw_3d,3); % time dimension

for iz =  
infile  = ['../data_after_manipulation/WvelTempGlobal_CESM_', ...
    num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];
% CESM wvel and temp 
w_ce_3d = ncread(infile, 'WVEL');
Tw_3d   = dummy.Tw_3d;
% ==================


%% === data analysis ===
for iday = 1 : ntime % loop through time
disp(['analyzing day ', num2str(iday,'%02d')])

% --- compute w', T', w'T' --- 
w_nl_prim = w_nl_3d(:,:,iday)-w_nl_av;
w_li_prim = w_li_3d(:,:,iday)-w_li_av;
w_ce_prim = w_ce_3d(:,:,iday)-w_ce_av;
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

% --- spatial av ---
[]=bootstrap5();
% ---

% --- clear vars to save memory --- 
clear w_nl_prim w_li_prim w_ce_prim Tw_prim 
clear wt_nl_prim wt_li_prim wt_ce_prim 
% ---------------------------------
end

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
