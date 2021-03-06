%% ===readme===

% descrip: matlab scripts plot Nonlinear Ekman vertical velocity (W)
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

  date_str='2020Feb20'; 

  header = 'This file was generated by DL via code ComputeSaveWvelTempGulfStreamDpthProfBorealSummerAv_2020Feb20.m';

infileW = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/WVEL/CESM_regrid_Apr-Sep_87-90_WVEL_Global.nc';
  ncdisp(infileW)
infileWav = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/WVEL/CESM_regrid_Apr-SepAv_WVEL_Global.nc';
  ncdisp(infileWav)

infileT = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/CESM_regrid_Apr-Sep_87-90_TEMP_Global.nc';
  ncdisp(infileT)
infileTav = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/CESM_regrid_Apr-SepAv_TEMP_Global.nc';
  ncdisp(infileTav)

outfileAv = ['../data_after_manipulation/TauKesaiRoWvelTempGlobalBorealSummerAv_',date_str,'.mat'];

% pic1 = ['../pics/TauKesaiRoEkmanWvelContoursGlobalBorealSummerAv_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))

% Constants4CESM_GulfStream

% tau_amp_csum  = 0;
  lon_1d_raw = [0.05   : 0.1 : 360];
  lat_1d_raw = [-89.95 : 0.1 : 90]';
  rho_w = 1020;
% --- Gulf Stream ---
% lat_limits = [-90 90]; % chosen lat ranges for comparison purpose
lat_limits = [25 50]; % chosen lat ranges
lon_limits = [275 325];
step_lat = 1;
step_lon = 1;
% vertical index, 9 elements totally, only use 50m depth for global pic
% same index for temp and wvel
start_z = 1;
count_z = 1;
stride_z = 1; 
% z_w_top = 0, 1000, 2000, 3000, 4000, 5000, 7000, 9000, 11000, [cm]
% z_t   = 500, 1500, 2500, 3500, 4500, 5500, 7500, 9500, 11500, [cm]
% ------------------
%=============================


%% === data loading and analysis ===
[~,indx4LatMin] = min(abs(lat_1d_raw-lat_limits(1)));% lat_1d_raw(indx4lat_limit1)
[~,indx4LatMax] = min(abs(lat_1d_raw-lat_limits(2)));% lat_1d_raw(indx4lat_limit2)
lat_1d = double(lat_1d_raw(indx4LatMin:step_lat:indx4LatMax));
count_lat = length(lat_1d);

[~,indx4LonMin] = min(abs(lon_1d_raw-lon_limits(1)));% lon_1d_raw(indx4lon_limit1)
[~,indx4LonMax] = min(abs(lon_1d_raw-lon_limits(2)));% lon_1d_raw(indx4lon_limit2)
lon_1d = double(lon_1d_raw(indx4LonMin:step_lon:indx4LonMax));
count_lon = length(lon_1d);


count_4dvar  = [count_lon   count_lat   count_z   Inf];
stride_4dvar = [step_lon    step_lat    stride_z  1];
start_4dvarAv  = [indx4LonMin indx4LatMin 1   1]; 

for iz = 1:9 % read one depth level one time
    
  disp(['working on ']) 
    
  start_4dvar  = [indx4LonMin indx4LatMin iz   1];   
 
  Tw_av = ncread(infileTav,'TEMP',start_4dvar,count_4dvar,stride_4dvar);
  Tw_av = ncread(infileTav,'TEMP',start_4dvar,count_4dvar,stride_4dvar);
 
wvel   = squeeze(ncread(infile_WVEL,'WVEL',  start_4dvar,count_4dvar,stride_4dvar)).*0.01;
z_wvel_ce = double(ncread(infile_WVEL,'z_w_top',start_z,count_z,stride_z));

infile_TEMP  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/' ...
    'CESM_regrid_MonthlyOutput_TEMP_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
disp(['loading ', infile_TEMP])
Tw   = squeeze(ncread(infile,'TEMP',start_4dvar,count_4dvar,stride_4dvar));
z_Tw = double(ncread(infile_TEMP,'z_t',  start_z,count_z,stride_z));


end

for iyr = 87:90
  for imon = [4:9]; % summer
      
%% === load data ===
infile_TTS  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TAUXTAUYSSH/' ...
    'CESM_regrid_MonthlyOutput_TAUXTAUYSSH_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
disp(['loading ', infile_TTS])
taux = ncread(infile_TTS,'TAUX', start_3dvar,count_3dvar,stride_3dvar).*0.1; 
tauy = ncread(infile_TTS,'TAUY', start_3dvar,count_3dvar,stride_3dvar).*0.1; 
ssh  = ncread(infile_TTS,'SSH',  start_3dvar,count_3dvar,stride_3dvar).*0.01;
time = ncread(infile_TTS,'time');
ntime= numel(time);

infile_WVEL  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/WVEL/' ...
    'CESM_regrid_MonthlyOutput_WVEL_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
disp(['loading ', infile_WVEL])
wvel   = squeeze(ncread(infile_WVEL,'WVEL',  start_4dvar,count_4dvar,stride_4dvar)).*0.01;
z_wvel_ce = double(ncread(infile_WVEL,'z_w_top',start_z,count_z,stride_z));

infile_TEMP  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TEMP/' ...
    'CESM_regrid_MonthlyOutput_TEMP_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
disp(['loading ', infile_TEMP])
Tw   = squeeze(ncread(infile,'TEMP',start_4dvar,count_4dvar,stride_4dvar));
z_Tw = double(ncread(infile_TEMP,'z_t',  start_z,count_z,stride_z));

outfile = ['../data_after_manipulation/' ...
    'WvelTempGlobal_CESM_',num2str(iyr),'-',num2str(imon,'%02d'),'.mat'];% ncdisp(infile)
% ==================


%% === data analysis ===
for iday = 1 : ntime % loop through time
  disp(['analyzing day ', num2str(iday,'%02d')])

% --- compute wind stress amplitude and curl ---
  % compute wind stress amplitude 
  tau_amp  = sqrt((taux(:,:,iday)').^2+(tauy(:,:,iday)').^2);   
  % compute wind stress curl
  tau_curl = CalcCurlz4UnevenGridsFunc(x_2d,y_2d, ...
     taux(:,:,iday)',tauy(:,:,iday)');
% ----------------------------------------------

% --- compute geostrophic velocity and geostrophic vorticity ---
  % compute geostrophic velocity
  [Ug,Vg] = CalcGeostrophyVel4UnevenGridsFunc( ...
       ssh(:,:,iday)',x_2d,y_2d,lat_1d); 
  % compute geostrophic vorticity
  kesai = CalcCurlz4UnevenGridsFunc(x_2d,y_2d, ...
           Ug,Vg);
       
  % compute Rossby num
  Ro = kesai./f_rp; 
% --------------------------------------------

% --- compute linear and nonlinear Ekman Wvel ---
% compute linear Ekman Wvel
  w_nl = CalcEkmanWvelFunc(rho_w,x_2d,y_2d, ...
     taux(:,:,iday)',tauy(:,:,iday)',f_rp,kesai,f_min);
% compute nonlinear Ekman Wvel 
  w_li = CalcEkmanWvelFunc(rho_w,x_2d,y_2d, ...
     taux(:,:,iday)',tauy(:,:,iday)',f_rp,0,f_min);
% df: difference between nonliner and linear Ekman Wvel
  w_df = w_nl - w_li; 
% 3d array for output files
  w_nl_3d(:,:,iday) = w_nl;
  w_li_3d(:,:,iday) = w_li;
  w_ce_3d(:,:,iday) = wvel(:,:,iday)';
  Tw_3d(:,:,iday)   = Tw(:,:,iday)';
% -------------------------

% --- compute cumulative sum (csum) --- 
  tau_amp_csum  = tau_amp_csum  + tau_amp;
  tau_curl_csum = tau_curl_csum + tau_curl;
  kesai_csum    = kesai_csum + kesai;
  Ro_csum       = Ro_csum + Ro;
  w_nl_csum     = w_nl_csum + w_nl;
  w_li_csum     = w_li_csum + w_li;
  w_df_csum     = w_df_csum + w_df;
  w_ce_csum     = w_ce_csum + wvel(:,:,iday)';
  Tw_csum       = Tw_csum   + Tw(:,:,iday)';
% ---------------

% --- compute number of loops ---
Ndays = Ndays + 1;
% -------------------------------

% --- clear vars to save memory --- 
clear tau_amp tau_curl Ug Vg kesai Ro w_nl w_li w_df 
% ---------------------------------
end

save(outfile,'header','lon_1d','lat_1d','w_nl_3d','w_li_3d', ...
    'w_ce_3d','z_wvel_ce','Tw_3d','z_Tw')

% --- clear vars to save memory --- 
  clear infile_TTS infile_WVEL infile_TEMP taux tauy ssh time ntime  
  clear wvel Tw w_nl_3d w_li_3d w_ce_3d Tw_3d
% ---------------------------------

  end
end

% --- compute mean value ---
tau_amp_TimeAv  = tau_amp_csum ./ Ndays;
tau_curl_TimeAv = tau_curl_csum ./ Ndays;
kesai_TimeAv    = kesai_csum ./ Ndays;
Ro_TimeAv       = Ro_csum ./ Ndays;
w_nl_TimeAv     = w_nl_csum ./ Ndays;
w_li_TimeAv     = w_li_csum ./ Ndays;
w_df_TimeAv     = w_df_csum ./ Ndays;
w_ce_TimeAv     = w_ce_csum ./ Ndays;
Tw_TimeAv       = Tw_csum ./ Ndays;
% --------------------------
% ==================
  
  
%% === output data ===
save(outfileAv,'header','lon_1d','lat_1d','Ndays','z_wvel_ce','z_Tw', ...
  'tau_amp_TimeAv','tau_curl_TimeAv','kesai_TimeAv','Ro_TimeAv', ...
  'w_nl_TimeAv','w_li_TimeAv','w_df_TimeAv','w_ce_TimeAv','Tw_TimeAv') 
% ====================
