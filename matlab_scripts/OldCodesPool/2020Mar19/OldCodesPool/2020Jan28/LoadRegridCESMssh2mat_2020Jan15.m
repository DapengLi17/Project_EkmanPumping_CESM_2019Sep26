%% ===readme===

% descrip: matlab scripts load ssh from CESM nc files, 
% convert the units, regrid and save data to .mat file   

% update history:
% v1.0 DL 2019Dec01

% extra notes:
% function required: LoadCESMKEncFiles3DVarFunc.m 
% =============


%% === set up environments ===
clear all;close all;clc;

date_str='2020Jan15';
addpath('Func4EkmanProject')
  
SSH_Str = {'SSH','time','TLONG','TLAT'};% SSH files
  
jultime = [datenum([0087,01,01]):datenum([0087,12,31]) ...
           datenum([0088,01,01]):datenum([0088,12,31]) ...
           datenum([0089,01,01]):datenum([0089,12,31]) ...
           datenum([0090,01,01]):datenum([0090,12,31])];
jultime(jultime==datenum([0088,02,29]))=[]; % delete 0088 (leap year) Feb 29

jultime_vec = datevec(jultime);

% --- input/output files ---
InDirSSH = '../raw_data/CESMncFilesKEFromAgartha_ZhangQiuying/SSH';
Outfile1 = ['../data_after_manipulation/CESMssh_KE_',date_str,'.mat'];
Pic1 = ['../pics/CheckRegrid/CheckRegridSSH_',date_str,'.png'];
% --------------------------

% --- lat/lon limits for regrid ---
% --- regrid data, interplate data from (u)t-grid to s-grid ---
% min(tlon(:)) = 129.1848; % t-129.18 --> s-130.0
% min(tlon(:)) = 129.2305

% max(tlon(:)) = 170.5412; % t-170.54 --> s-170.0
% max(tlon(:)) = 170.5952  

% min(tlat(:)) = 25.5623;% t- 25.56 --> s- 26.0
% min(tlat(:)) = 25.6074 

% max(tlat(:)) = 46.1919;% t- 46.19 --> s- 46.0
% max(tlat(:)) = 46.2335
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
% -------------------------------
% =======================


%% === load data ===
[~,tlat,tlon,ssh_raw] = LoadCESMKEncFiles3DVarFunc(InDirSSH,SSH_Str);
% ===================


%% === data analysis ===
% unit of ssh is centimeter

ssh_unitconvt = ssh_raw./100; % [m]

% % check raw nc files for comparison
% dummy_infile = '../raw_data/CESMncFilesKEFromAgartha_ZhangQiuying/SSH/kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_SSH_00870101-00871231.nc';
% ncdisp(dummy_infile)
% Variables:
%     SSH       
%            Size:       401x226x365
%            Dimensions: nlon,nlat,time
%            Datatype:   single
%            Attributes:
%                        long_name     = 'Sea Surface Height'
%                        units         = 'centimeter'
%                        coordinates   = 'TLONG TLAT time'
%                        grid_loc      = '2110'
%                        cell_methods  = 'time: mean'
%                        _FillValue    = 9.969209968386869e+36
%                        missing_value = 9.969209968386869e+36


% --- regrid data, interplate data from (u)t-grid to s-grid ---
tlat_rot = rot90(tlat);
tlon_rot = rot90(tlon);

for iDay = 1:size(ssh_unitconvt,3) % loop through time
   ssh_rot(:,:,iDay) = rot90(ssh_unitconvt(:,:,iDay));
end

[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

for iDay=1:size(ssh_rot,3) % loop through time
  iDay
  ssh(:,:,iDay) = griddata(tlon_rot,tlat_rot, ...
      ssh_rot(:,:,iDay),Lon,Lat);
end
% ====================


%% === make pics ===
f1=figure; % check rotation

  set(f1,'units','normalized','position',[0,0,1,1]);
  plt_indx = 374;
  ssh_plt = ssh_unitconvt(:,:,plt_indx);
  ssh_pltmin=nanmin(ssh_plt(:));
  ssh_pltmax=nanmax(ssh_plt(:));  
  
 subplot(2,2,1) % raw data
  pcolor(tlon,tlat,ssh_unitconvt(:,:,plt_indx));shading interp;
  caxis([ssh_pltmin ssh_pltmax]);hc=colorbar;title(hc,'[N/m2]')
  title(['raw CESM ssh, time: ',datestr(jultime(plt_indx))]);
  
 subplot(2,2,2) % data after rotation
  pcolor(tlon_rot,tlat_rot,ssh_rot(:,:,plt_indx));shading interp;
  caxis([ssh_pltmin ssh_pltmax]);hc=colorbar;title(hc,'[N/m2]')
  title(['rotated CESM ssh, time: ',datestr(jultime(plt_indx))]);
  
 subplot(2,2,3) % data after rotation and interpolation
  pcolor(Lon,Lat,ssh(:,:,plt_indx));shading interp;
  caxis([ssh_pltmin ssh_pltmax]);hc=colorbar;title(hc,'[N/m2]')
  title(['rotated and interpolated CESM ssh,  time: ', ...
      datestr(jultime(plt_indx))]);
 
 subplot(2,2,4) % time series
  ssh_llav = squeeze(nanmean(nanmean(abs(ssh),1),2));
  ssh_llmax = squeeze(nanmax(nanmax(abs(ssh),[],1),[],2));  
  h1=plot(jultime,ssh_llav,'b');hold on;
  h2=plot(jultime,ssh_llmax,'r');grid on;
    legend([h1,h2],{'lat-lon av ssh','lat-lon max ssh'});
    set(gca,'ylim',[0 1]);ylabel('N/m2');datetick('x','mmm/yy');
    title(['rotated and interpolated CESM ssh time series', ...
        datestr(jultime(plt_indx))]);  
% ==============


%% === save data ===
print(f1,'-dpng',Pic1)

header = ['This file was created by DL on ',date_str,' via code ', ...
    'ConvertCESMssh2mat_2020Jan15.m, units of ssh is ', ...
    'converted to m.'];
save(Outfile1,'header','jultime','jultime_vec','tlon','tlat', ...
    'ssh','-v7');
%===================
