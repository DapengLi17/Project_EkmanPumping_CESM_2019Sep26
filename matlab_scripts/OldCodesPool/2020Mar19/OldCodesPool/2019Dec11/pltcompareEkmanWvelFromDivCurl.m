%% ===readme===

% descrip: matlab scripts plot tau, UVW vel, Kesai from 
% 1) CESM model output
% 2) Stern nonlinear Ekman vertical pumping velocity
% 3) linear Ekman vertical pumping velocity

% update history:
% v1.0 DL 2019Nov14

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Nov14';
  
infile1 = '../data_after_manipulation/CESMtaux_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMtauy_KE_2019Dec01.mat';  
infile3 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
infile4 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';

addpath(genpath('Func4EkmanProject/'))
addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))
rho0 = 1020;
%=============================


%% === load data ===
TAUX_CE = load(infile1);
TAUY_CE = load(infile2);
UVEL_CE = load(infile3);
VVEL_CE = load(infile4);
%===================


%% === data analysis ===
% convert lat and lon to x and y
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

% Lat_vec=U_CE.Lat_vec;Lon_vec=U_CE.Lon_vec;
% Lon=U_CE.Lon;Lat=U_CE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);

uvel_raw = UVEL_CE.u;
vvel_raw = VVEL_CE.v;
% w = U_CE.w;
taux_raw = TAUX_CE.taux;
tauy_raw = TAUY_CE.tauy;
% tau_A = sqrt(taux.^2+tauy.^2);
% tauA_av = squeeze(nanmean(nanmean(tauA,1),2));
jultime_vec = UVEL_CE.jultime_vec;
jultime= datenum(jultime_vec);

clear TAUX_CE TAUY_CE UVEL_CE VVEL_CE

f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
for i = 1 : length(jultime) % loop through time
    kesai(:,:,i) = CalcRelVorticityFunc(x,y,uvel_raw(:,:,4,i),vvel_raw(:,:,4,i));
end

% compute nonlinear Ekman Wvel from curl and divergence Eq
for i = 1 : length(jultime) % loop through time
 [w_cur(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho0,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,kesai(:,:,i));
 [w_div(:,:,i)] = CalcEkmanWvelStern65Func(rho0,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,kesai(:,:,i));
end


%%
stime = datenum([0088 01 01 00 00 00]);
ftime = datenum([0088 01 01 23 59 59]);
IndxTime = find(jultime>=stime & jultime<=ftime);

%%
figure;
 
 w_pltmin = -1e-4;
 w_pltmax = 1e-4;
 
 subplot(2,2,1)
  pcolor(Lon,Lat,w_cur(:,:,51));shading interp;
%   [c,h]=contourf(Lon,Lat,w_cur(:,:,51));set(h,'linestyle','none');
  caxis([w_pltmin w_pltmax]);polarmap(64,1);colorbar; 
%   set(gca,'xlim',[120 169.5],'ylim',[26.5 55]);
  title(['Rossby Num: kesai/f']);
 subplot(2,2,2)
  pcolor(Lon,Lat,w_div(:,:,51));shading interp;
%   [c,h]=contourf(Lon,Lat,Ro_SeaStd(:,:,1),linspace(-0.4,0.4,20));set(h,'linestyle','none');
  caxis([w_pltmin w_pltmax]);polarmap(64,1);colorbar; 
%   set(gca,'xlim',[120 169.5],'ylim',[26.5 55]);
  title(['Rossby Num: kesai/f']);
 subplot(2,2,3)
  pcolor(Lon,Lat,w_cur(:,:,366));shading interp;
%   [c,h]=contourf(Lon,Lat,w_cur(:,:,51));set(h,'linestyle','none');
  caxis([w_pltmin w_pltmax]);polarmap(64,1);colorbar; 
%   set(gca,'xlim',[120 169.5],'ylim',[26.5 55]);
  title(['Rossby Num: kesai/f']);
 subplot(2,2,4)
  pcolor(Lon,Lat,w_div(:,:,366));shading interp;
%   [c,h]=contourf(Lon,Lat,Ro_SeaStd(:,:,1),linspace(-0.4,0.4,20));set(h,'linestyle','none');
  caxis([w_pltmin w_pltmax]);polarmap(64,1);colorbar; 
%   set(gca,'xlim',[120 169.5],'ylim',[26.5 55]);
  title(['Rossby Num: kesai/f']);
