%% ===readme===

% descrip: matlab scripts load u from CESM nc files, 
% convert the units, regrid and save data to .mat file   

% update history:
% v1.0 DL 2019Dec01

% extra notes:
% function required: LoadCESMKEncFiles3DVarFunc.m 
% =============


%% === set up environments ===
clear all;close all;clc;

date_str='2019Dec01';
addpath('Func4EkmanProject')
  
UVEL_Str = {'UVEL','time','ULONG','ULAT','z_t'}; % UVEL files
  
jultime = [datenum([0087,01,01]):datenum([0087,12,31]) ...
           datenum([0088,01,01]):datenum([0088,12,31]) ...
           datenum([0089,01,01]):datenum([0089,12,31]) ...
           datenum([0090,01,01]):datenum([0090,12,31])];
jultime(jultime==datenum([0088,02,29]))=[]; % delete 0088 (leap year) Feb 29

jultime_vec = datevec(jultime);

% --- input/output files ---
InDirUVEL = '../raw_data/CESMncFilesKEFromAgartha_ZhangQiuying/UVEL';
Outfile1 = ['../data_after_manipulation/CESMuvel_KE_',date_str,'.mat'];
Pic1 = ['../pics/CheckRegrid/CheckRegridUVEL_',date_str,'.png'];
% --------------------------

% --- lat/lon limits for regrid ---
% --- regrid data, interplate data from (u)t-grid to s-grid ---
% min(tlon(:)) = 129.1848; % t-129.18 --> s-130.0
% min(ulon(:)) = 129.2305

% max(tlon(:)) = 170.5412; % t-170.54 --> s-170.0
% max(ulon(:)) = 170.5952  

% min(tlat(:)) = 25.5623;% t- 25.56 --> s- 26.0
% min(ulat(:)) = 25.6074 

% max(tlat(:)) = 46.1919;% t- 46.19 --> s- 46.0
% max(ulat(:)) = 46.2335
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
% -------------------------------
% =======================


%% === load data ===
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

IndxZ = 3 : 7;
  [~,ulat,ulon,u_raw,z_raw] = LoadCESMKEncFiles4DVarFunc( ...
        InDirUVEL,UVEL_Str,IndxZ);
% ===================


%% === data analysis ===
% unit of wind stress u and u are dyne/centimeter??
% 1 dyne/centimeter?? [dyn/cm??] = 0.1 newton/meter?? [N/m??]
% see https://www.translatorscafe.com/unit-converter/en/pressure/25-18/dyne%2Fcentimeter%C2%B2-newton%2Fmeter%C2%B2/

% cm/s to m/s
u_unitconvt = u_raw./100;
% centimeters to meters
z = double(z_raw./100);

% % check raw nc files for comparison
% dummy_infile = '../raw_data/CESMncFilesKEFromAgartha_ZhangQiuying/UVEL/kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_UVEL_00870101-00871231.nc';
% ncdisp(dummy_infile)
% Variables:
%     UVEL      
%            Size:       401x226x21x365
%            Dimensions: nlon,nlat,z_t,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Velocity in grid-x direction'
%                        units         = 'centimeter/s'
%                        coordinates   = 'ULONG ULAT z_t time'
%                        grid_loc      = '3221'
%                        cell_methods  = 'time: mean'
%                        _FillValue    = 9.969209968386869e+36
%                        missing_value = 9.969209968386869e+36

% --- regrid data, interplate data from (u)t-grid to s-grid ---
ulat_rot = rot90(ulat);
ulon_rot = rot90(ulon);
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

for iDay = 1 : size(u_unitconvt,4) % loop through time
    iDay
    for iZ = 1 : size(u_unitconvt,3) % loop through depth
      u_rot(:,:,iZ,iDay) = rot90(u_unitconvt(:,:,iZ,iDay));
      u(:,:,iZ,iDay) = griddata(ulon_rot,ulat_rot, ...
         u_rot(:,:,iZ,iDay),Lon,Lat);
    end
end
% ====================


%% === make pics ===

% ncread(dummy_infile,'UVEL')

f1=figure; % check rotation

  set(f1,'units','normalized','position',[0,0,1,1]);
  pltTimeIndx = 51;
  pltZIndx = 3;
  u_plt = u_unitconvt(:,:,pltZIndx,pltTimeIndx);
  u_pltmin=nanmin(u_plt(:));
  u_pltmax=nanmax(u_plt(:));  

 subplot(2,2,1) %
  u_llav = squeeze(nanmean(nanmean(abs(u),1),2));
  u_llmax = squeeze(nanmax(nanmax(abs(u),[],1),[],2));  
  h1=plot(jultime,u_llav,'b');hold on;
  h2=plot(jultime,u_llmax,'r');grid on;
    legend([h1,h2],{'lat-lon av u','lat-lon max u'});
    set(gca,'ylim',[0 1]);ylabel('N/m2');datetick('x','mmm/yy');
    title(['rotated and interpolated CESM u time series', ...
        datestr(jultime(pltTimeIndx))]);   
  
 subplot(2,2,2) % raw data
  pcolor(ulon,ulat,u_unitconvt(:,:,pltZIndx,pltTimeIndx));shading interp;
  caxis([u_pltmin u_pltmax]);hc=colorbar;title(hc,'[m/s]')
  title(['raw CESM u, time: ',datestr(jultime(pltTimeIndx))]);
  
 subplot(2,2,3) % data after rotation
  pcolor(ulon_rot,ulat_rot,u_rot(:,:,pltZIndx,pltTimeIndx));shading interp;
  caxis([u_pltmin u_pltmax]);hc=colorbar;title(hc,'[m/s]')
  title(['rotated CESM u, time: ',datestr(jultime(pltTimeIndx))]);
  
 subplot(2,2,4) % data after rotation and interpolation
  pcolor(Lon,Lat,u(:,:,pltZIndx,pltTimeIndx));shading interp;
  caxis([u_pltmin u_pltmax]);hc=colorbar;title(hc,'[m/s]')
  title(['rotated and interpolated CESM u,  time: ', ...
      datestr(jultime(pltTimeIndx))]);
 
% ==============


%% === save data ===
print(f1,'-dpng',Pic1)

header = ['This file was created by DL on ',date_str,' via code ', ...
    'ConvertCESMuvel2mat_2019Dec01.m, units of u is ', ...
    'converted to m/s.'];
save(Outfile1,'header','jultime','jultime_vec','ulon','ulat', ...
    'u','z','-v7.3');
%===================
