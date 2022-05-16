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
  
infile1 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';
pic1    = ['../pics/RossbyNum_KE_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))
addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))
rho0 = 1020;
%=============================


%% === load data ===
U_CE = load(infile1);
V_CE = load(infile2);
%===================


%% === data analysis ===
% convert lat and lon to x and y
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

% Lat_vec=U_CE.Lat_vec;Lon_vec=U_CE.Lon_vec;
% Lon=U_CE.Lon;Lat=U_CE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);

u_raw = U_CE.u;
v_raw = V_CE.v;
% w = U_CE.w;
% taux = U_CE.taux;
% tauy = U_CE.tauy;
% tau_A = sqrt(taux.^2+tauy.^2);
% tauA_av = squeeze(nanmean(nanmean(tauA,1),2));
jultime_vec = U_CE.jultime_vec;
jultime= datenum(jultime_vec);
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
U_CE.z(4)
for i = 1 : length(jultime) % loop through time
    kesai(:,:,i) = CalcRelVorticityFunc(x,y,u_raw(:,:,4,i),v_raw(:,:,4,i));
end

% compute Rossby number
for i = 1 : length(jultime)
    Ro(:,:,i) = kesai(:,:,i)./f;  
end

stime = datenum([0088 01 01 00 00 00]);
ftime = datenum([0088 01 01 23 59 59]);
IndxTime = find(jultime>=stime & jultime<=ftime);
 
IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);

% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{12}]; % Winter
  IndxSea{2} = [IndxMon{3};IndxMon{4};IndxMon{5}]; % Spring 
  IndxSea{3} = [IndxMon{6};IndxMon{7};IndxMon{8}]; % Summer
  IndxSea{4} = [IndxMon{9};IndxMon{10};IndxMon{11}]; % Fall
  
Ro_Sea = Ro(:,:,[IndxMon{1};IndxMon{2};IndxMon{3}; ...
    IndxMon{10};IndxMon{11};IndxMon{12}]);
% Ro_Sea = Ro(:,:,[335:424]);
Ro_SeaStd = nanstd(Ro_Sea,0,3);

%% === make pics ===
f1=figure;
 subplot(1,2,1)
%   pcolor(Lon,Lat,Ro(:,:,IndxTime));shading interp;
  [c,h]=contourf(Lon,Lat,Ro(:,:,IndxTime),linspace(-0.4,0.4,20));set(h,'linestyle','none');
  caxis([-0.4 0.4]);polarmap(64,1);colorbar; 
  set(gca,'xlim',[120 200],'xtick',[120:15:210], ...
      'ylim',[24 60],'ytick',[24:6:60]);
  title(['Snapshot of Rossby Num (\xi/f) on 0088Jan01']);
  
 subplot(1,2,2)
%   pcolor(Lon,Lat,Ro_SeaStd(:,:,1));shading interp;
  [c,h]=contourf(Lon,Lat,Ro_SeaStd(:,:,1),linspace(-0.4,0.4,20));set(h,'linestyle','none');
  polarmap;colorbar;caxis([0 0.4]);
  set(gca,'xlim',[120 200],'xtick',[120:15:210], ...
      'ylim',[24 60],'ytick',[24:6:60]);
  title(['RMS of Rossby Num (\xi/f) in winter seasons']);
% ==================


%% === output data ===
print(f1,'-dpng',pic1)
% ====================