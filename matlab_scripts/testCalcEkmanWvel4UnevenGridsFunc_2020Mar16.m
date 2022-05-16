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

  date_str='2020Mar16'; 

pic1 = ['../pics/testCalcEkmanWvel4UnevenGridsFunc_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))

Constants4CESM_Global
%=============================

      
%% === load data ===
iyr = 88;
imon = 01; 
infile_TTS  = ['../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TAUXTAUYSSH/' ...
    'CESM_regrid_MonthlyOutput_TAUXTAUYSSH_',num2str(iyr),'-',num2str(imon,'%02d'),'.nc'];% ncdisp(infile)
disp(['loading ', infile_TTS])
taux = ncread(infile_TTS,'TAUX', start_3dvar,count_3dvar,stride_3dvar).*0.1; 
tauy = ncread(infile_TTS,'TAUY', start_3dvar,count_3dvar,stride_3dvar).*0.1; 
ssh  = ncread(infile_TTS,'SSH',  start_3dvar,count_3dvar,stride_3dvar).*0.01;
time = ncread(infile_TTS,'time');
ntime= numel(time);
% ==================


%% === data analysis ===
for iday = 1 : ntime % loop through time
  disp(['analyzing day ', num2str(iday,'%02d')])

% --- compute wind stress amplitude and curl ---
%   % compute wind stress amplitude 
%   tau_amp  = sqrt((taux(:,:,iday)').^2+(tauy(:,:,iday)').^2);   
%   % compute wind stress curl
%   tau_curl = CalcCurlz4UnevenGridsFunc(x_2d,y_2d, ...
%      taux(:,:,iday)',tauy(:,:,iday)');
% ----------------------------------------------

% --- compute geostrophic velocity and geostrophic vorticity ---
  % compute geostrophic velocity
  [Ug,Vg] = CalcGeostrophyVel4UnevenGridsFunc( ...
       ssh(:,:,iday)',x_2d,y_2d,lat_1d); 
  % compute geostrophic vorticity
  kesai = CalcCurlz4UnevenGridsFunc(x_2d,y_2d, ...
           Ug,Vg);
       
  % compute Rossby num
  Ro(:,:,iday) = kesai./f_rp; 
  
  % --- compute linear and nonlinear Ekman Wvel ---
  [W_TE(:,:,iday), W_LE(:,:,iday), W_NE(:,:,iday)] = ...
      CalcEkmanWvel4UnevenGridsFunc(rho_w, x_2d, y_2d, ...
      taux(:,:,iday)', tauy(:,:,iday)', f_rp, kesai);
  
end
% ==================


%% === make pics ===  
f1=figure;
 
 iday = 5;
 set(f1,'units','normalized','position',[0,0,1,1])
 m_proj('miller','lon',lon_limits,'lat',lat_limits);

 subplot(2,2,1)
   m_pcolor(lon_1d,lat_1d,W_LE(:,:,iday).*86400);shading interp;
    polarmap;caxis([-2 2]);colorbar;
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none');axis normal;
   title('Linear Ekman Pumping [m/day]');

 subplot(2,2,2)
   m_pcolor(lon_1d,lat_1d,W_NE(:,:,iday).*86400);shading interp;
    polarmap;caxis([-2 2]);colorbar;
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none');axis normal;
   title('Non-Linear Ekman Pumping [m/day]');
   
 subplot(2,2,3)
   W_sum = W_LE+W_NE;
   m_pcolor(lon_1d,lat_1d,W_sum(:,:,iday).*86400);shading interp;
    polarmap;caxis([-2 2]);colorbar;
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none');axis normal; 
   title('Linear + Non-Linear Ekman Pumping [m/day]');
   
 subplot(2,2,4)
   m_pcolor(lon_1d,lat_1d,W_TE(:,:,iday).*86400);shading interp;
    polarmap;caxis([-2 2]);colorbar;
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none');axis normal;
   title('Total Ekman Pumping [m/day]');

% --- plot Rossub Num ---
% Ro_plt = Ro(:,:,1);
% Ro_plt = mean(Ro,3);
% indx_Lat= find(abs(lat_1d)<10);
% Ro_plt(indx_Lat,:) = nan;

% f1=figure;
%   m_proj('miller','lon',lon_limits,'lat',lat_limits);
%   m_pcolor(lon_1d,lat_1d,Ro_plt);shading interp;
%   polarmap; caxis([-1 1]); colorbar;
%   m_coast('patch',[.7 .7 .7],'edgecolor','none');
%   m_grid('linestyle','none');axis normal;
%   title('Rossby Num. 88-01 av');
%   print(gcf,'-dpng','test.png');
% ----------------------
% =======================
   
   
%% === save pics ===
print(f1,'-dpng',pic1);
% ==================