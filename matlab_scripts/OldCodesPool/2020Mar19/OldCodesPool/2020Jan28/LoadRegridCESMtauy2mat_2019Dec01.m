%% ===readme===

% descrip: matlab scripts load tauy from CESM nc files, 
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
  
TAUY_Str = {'TAUY','time','ULONG','ULAT'};% TAUX files
  
jultime = [datenum([0087,01,01]):datenum([0087,12,31]) ...
           datenum([0088,01,01]):datenum([0088,12,31]) ...
           datenum([0089,01,01]):datenum([0089,12,31]) ...
           datenum([0090,01,01]):datenum([0090,12,31])];
jultime(jultime==datenum([0088,02,29]))=[]; % delete 0088 (leap year) Feb 29

jultime_vec = datevec(jultime);

% --- input/output files ---
InDirTAUY = '../raw_data/CESMncFilesKEFromAgartha_ZhangQiuying/TAUY';
Outfile1 = ['../data_after_manipulation/CESMtauy_KE_',date_str,'.mat'];
Pic1 = ['../pics/CheckRegrid/CheckRegridTAUY_',date_str,'.png'];
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
[~,ulat,ulon,tauy_raw] = LoadCESMKEncFiles3DVarFunc(InDirTAUY,TAUY_Str);
% ===================


%% === data analysis ===
% unit of wind stress tauy and tauy are dyne/centimeter²
% 1 dyne/centimeter² [dyn/cm²] = 0.1 newton/meter² [N/m²]
% see https://www.translatorscafe.com/unit-converter/en/pressure/25-18/dyne%2Fcentimeter%C2%B2-newton%2Fmeter%C2%B2/

tauy_unitconvt = tauy_raw.*0.1; % [N/m²]

% % check raw nc files for comparison
% dummy_infile = '../raw_data/CESMncFilesKEFromAgartha_ZhangQiuying/TAUX/kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_TAUX_00870101-00871231.nc';
% ncdisp(dummy_infile)
% Variables:
%     TAUX      
%            Size:       401x226x365
%            Dimensions: nlon,nlat,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Windstress in grid-x direction'
%                        units         = 'dyne/centimeter^2'
%                        coordinates   = 'ULONG ULAT time'
%                        grid_loc      = '2220'
%                        cell_methods  = 'time: mean'
%                        _FillValue    = 9.969209968386869e+36
%                        missing_value = 9.969209968386869e+36


% --- regrid data, interplate data from (u)t-grid to s-grid ---
ulat_rot = rot90(ulat);
ulon_rot = rot90(ulon);

for iDay = 1:size(tauy_unitconvt,3) % loop through time
   tauy_rot(:,:,iDay) = rot90(tauy_unitconvt(:,:,iDay));
end

[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

for iDay=1:size(tauy_rot,3) % loop through time
  iDay
  tauy(:,:,iDay) = griddata(ulon_rot,ulat_rot, ...
      tauy_rot(:,:,iDay),Lon,Lat);
end
% ====================


%% === make pics ===
f1=figure; % check rotation

  set(f1,'units','normalized','position',[0,0,1,1]);
  plt_indx = 374;
  tauy_plt = tauy_unitconvt(:,:,plt_indx);
  tauy_pltmin=nanmin(tauy_plt(:));
  tauy_pltmax=nanmax(tauy_plt(:));  
  
 subplot(2,2,1) % raw data
  pcolor(ulon,ulat,tauy_unitconvt(:,:,plt_indx));shading interp;
  caxis([tauy_pltmin tauy_pltmax]);hc=colorbar;title(hc,'[N/m2]')
  title(['raw CESM tauy, time: ',datestr(jultime(plt_indx))]);
  
 subplot(2,2,2) % data after rotation
  pcolor(ulon_rot,ulat_rot,tauy_rot(:,:,plt_indx));shading interp;
  caxis([tauy_pltmin tauy_pltmax]);hc=colorbar;title(hc,'[N/m2]')
  title(['rotated CESM tauy, time: ',datestr(jultime(plt_indx))]);
  
 subplot(2,2,3) % data after rotation and interpolation
  pcolor(Lon,Lat,tauy(:,:,plt_indx));shading interp;
  caxis([tauy_pltmin tauy_pltmax]);hc=colorbar;title(hc,'[N/m2]')
  title(['rotated and interpolated CESM tauy,  time: ', ...
      datestr(jultime(plt_indx))]);
 
 subplot(2,2,4) % time series
  tauy_llav = squeeze(nanmean(nanmean(abs(tauy),1),2));
  tauy_llmax = squeeze(nanmax(nanmax(abs(tauy),[],1),[],2));  
  h1=plot(jultime,tauy_llav,'b');hold on;
  h2=plot(jultime,tauy_llmax,'r');grid on;
    legend([h1,h2],{'lat-lon av tauy','lat-lon max tauy'});
    set(gca,'ylim',[0 1]);ylabel('N/m2');datetick('x','mmm/yy');
    title(['rotated and interpolated CESM tauy time series', ...
        datestr(jultime(plt_indx))]);  
% ==============


%% === save data ===
print(f1,'-dpng',Pic1)

header = ['This file was created by DL on ',date_str,' via code ', ...
    'ConvertCESMtauy2mat_2019Dec01.m, units of tauy is ', ...
    'converted to N/m2.'];
save(Outfile1,'header','jultime','jultime_vec','ulon','ulat', ...
    'tauy','-v7');
%===================
