%% ===readme===

% descrip: matlab scripts extract parameters from CESM nc files, 
% regrid the data and save them to .mat file   

% update history:
% v1.0 DL 2019Oct05

% extra notes:
% parameters needed: u,v,w velocity (UVEL,VVEL,WVEL), 
% wind stress (TAUX, TAUY) and temperature (TEMP). 
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Nov22';
  addpath('Func4EkmanProject')
  
  TAUX_Str = {'TAUX','time','ULONG','ULAT'};% TAUX files
  TAUY_Str = {'TAUY','time','ULONG','ULAT'}; % TAUY files
  UVEL_Str = {'UVEL','time','ULONG','ULAT','z_t'}; % UVEL files
  VVEL_Str = {'VVEL','time','ULONG','ULAT','z_t'}; % VVEL files
  WVEL_Str = {'WVEL','time','TLONG','TLAT','z_w_top'}; % WVEL files
  TEMP_Str = {'TEMP','time','TLONG','TLAT','z_t'};
  
jultime = datenum([0087,01,01]):datenum([0087,12,31]);
jultime_vec = datevec(jultime);

IndxZu = 3; 
IndxZw = 3;
% check depth of u,v,w, temp
% ncdump -v z_w_top kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_WVEL_00870101-00871231.nc 
% unit: cm
%  z_w_top = 0, 1000, 2000, 3000, 4000, 5000, 7000, 9000, 11000, 13000, 15000, 
%     19182.12, 24358.45, 29481.47, 33831.23, 39258.05, 50370.69, 60699.67, 
%     74439.8, 92804.35, 117104 ;

% ncdump -v z_t kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_UVEL_00870101-00871231.nc 
% ncdump -v z_t kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_TEMP_00870101-00871231.nc 
% unit: cm
%  z_t = 500, 1500, 2500, 3500, 4500, 5500, 7500, 9500, 11500, 13500, 15500, 
%     19766.03, 25137.02, 30511.92, 35109.35, 40878.46, 52772.8, 63886.26, 
%     78700.25, 98470.59, 124456.7 ; 

%InDir = '/atlantis3/zhangqiuying/regrid_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02/ocn_kuroshio/';
InDir1 = '../raw_data/';
Outfile1 = ['../data_after_manipulation/CESM_TauUVWTw_dpthindxof', ...
    num2str(IndxZu),'_KE_',date_str,'.mat'];
Pic1 = ['../pics/Check2DInterpolateUVel_',date_str,'.png'];
%=============================


%% === load data ===
[~,ulat,ulon,taux_raw] = LoadCESMKEncFiles3DVarFunc(InDir1,TAUX_Str);
[~,~,~,tauy_raw] = LoadCESMKEncFiles3DVarFunc(InDir1,TAUY_Str);
[~,~,~,u_raw,z_u] = LoadCESMKEncFiles4DVarFunc(InDir1,UVEL_Str,IndxZu);
[~,~,~,v_raw,z_v] = LoadCESMKEncFiles4DVarFunc(InDir1,VVEL_Str,IndxZu);
[~,tlat,tlon,w_raw,z_w] = LoadCESMKEncFiles4DVarFunc(InDir1,WVEL_Str,IndxZw);
[~,~,~,temp_raw,z_t] = LoadCESMKEncFiles4DVarFunc(InDir1,TEMP_Str,IndxZu);
%===================


%% === data analysis ===
% --- regrid data, interplate data from (u)t-grid to s-grid ---
% min(tlon(:)) = 129.1848; % t-129.18 --> s-130.0
% min(ulon(:)) = 129.2305

% max(tlon(:)) = 170.5412; % t-170.54 --> s-170.0
% max(ulon(:)) = 170.5952  

% min(tlat(:)) = 25.5623;% t- 25.56 --> s- 26.0
% min(ulat(:)) = 25.6074 

% max(tlat(:)) = 46.1919;% t- 46.19 --> s- 46.0
% max(ulat(:)) = 46.2335
ulat_rot = rot90(ulat);
ulon_rot = rot90(ulon);
tlat_rot = rot90(tlat);
tlon_rot = rot90(tlon);

% centimeters to meters
dpth_u = double(z_u./100);
dpth_v = double(z_v./100);
dpth_w = double(z_w./100); 
dpth_temp = double(z_t./100); 

% unit of wind stress taux and tauy are dyne/centimeter²
% 1 dyne/centimeter² [dyn/cm²] = 0.1 newton/meter² [N/m²]
% see https://www.translatorscafe.com/unit-converter/en/pressure/25-18/dyne%2Fcentimeter%C2%B2-newton%2Fmeter%C2%B2/
taux_unitconvt = taux_raw.*0.1;
tauy_unitconvt = tauy_raw.*0.1;
% cm/s to m/s
u_unitconvt = u_raw./100;
v_unitconvt = v_raw./100;
w_unitconvt = w_raw./100;

for iDay=1:size(taux_unitconvt,3)
   taux_rot(:,:,iDay) = rot90(taux_unitconvt(:,:,iDay));
   tauy_rot(:,:,iDay) = rot90(tauy_unitconvt(:,:,iDay));
   u_rot(:,:,iDay) = rot90(u_unitconvt(:,:,iDay)); 
   v_rot(:,:,iDay) = rot90(v_unitconvt(:,:,iDay));
   w_rot(:,:,iDay) = rot90(w_unitconvt(:,:,iDay));
   temp_rot(:,:,iDay) = rot90(temp_raw(:,:,iDay));
end

Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

for iDay=1:size(taux_rot,3)
  iDay
  taux(:,:,iDay) = griddata(ulon_rot,ulat_rot, ...
      taux_rot(:,:,iDay),Lon,Lat);
  tauy(:,:,iDay) = griddata(ulon_rot,ulat_rot, ...
      tauy_rot(:,:,iDay),Lon,Lat);
  u(:,:,iDay) = griddata(ulon_rot,ulat_rot,u_rot(:,:,iDay),Lon,Lat);
  v(:,:,iDay) = griddata(ulon_rot,ulat_rot,v_rot(:,:,iDay),Lon,Lat);
  w(:,:,iDay) = griddata(tlon_rot,tlat_rot,w_rot(:,:,iDay),Lon,Lat);
  temp(:,:,iDay) = griddata(tlon_rot,tlat_rot, ...
      temp_rot(:,:,iDay),Lon,Lat);
end
%=======================
 

%% === make pics ===
f1=figure; % check rotation
  set(f1,'units','normalized','position',[0,0.3,1,0.5]);
  plt_indx = 10;
 subplot(1,3,1) % raw data
  [c,h] = contourf(ulon,ulat,u_unitconvt(:,:,plt_indx));
  caxis([-1 2]);title('raw CESM UVel data');
 subplot(1,3,2) % data after rotation
  [c,h] = contourf(ulon_rot,ulat_rot,u_rot(:,:,plt_indx));
  caxis([-1 2]);title('rotated CESM UVel data');
 subplot(1,3,3) % data after rotation and interpolation
  [c,h] = contourf(Lon,Lat,u(:,:,plt_indx));
  caxis([-1 2]);colorbar;title('rotated and interpolated CESM UVel data');
% ==============


%% === save data ===
print(f1,'-dpng',Pic1)

header = ['This file was created by DL on ',date_str,' via code ', ...
    'ConvertCESMnc2mat_',date_str,'.m.'];
save(Outfile1,'header','jultime','jultime_vec','Lon','Lat','Lon_vec', ...
    'Lat_vec','taux','tauy','u','v','w','temp','dpth_u','dpth_v', ...
    'dpth_w','dpth_temp','-v7');
%===================
