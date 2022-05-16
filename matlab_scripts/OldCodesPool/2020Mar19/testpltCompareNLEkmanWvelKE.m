%% ===readme===

% descrip: matlab scripts plot composite-av W vel during winter storms from 
% 2) Stern nonlinear Ekman vertical pumping velocity
% 1) linear Ekman vertical pumping velocity
% 3) CESM model output

% update history:
% v1.0 DL 2019Dec05
% v1.1 DL 2019Dec09

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2020Jan15';
  
infile1 = '../data_after_manipulation/CESMtaux_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMtauy_KE_2019Dec01.mat';  
% infile3 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
% infile4 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';
infile5 = '../data_after_manipulation/CESMwvel_KE_2019Dec01.mat';
infile7 = '../data_after_manipulation/CESMssh_KE_2020Jan15.mat';

pic1 = ['../pics/SpaceAvWindStress_KE_',date_str,'.png'];
pic2 = ['../pics/WvelContours_KE_',date_str,'.png'];
pic3 = ['../pics/Coastline_KE_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))
addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))
rho_w = 1020;

% lat and lon range
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
Ro_std_limit = 0.15; 
% tau_st_limt = 0.25; 
tau_st_limt = 0.2;
indxZ = 3;
%=============================


%% === load data ===
TAUX_CE = load(infile1);
TAUY_CE = load(infile2);
% UVEL_CE = load(infile3);
% VVEL_CE = load(infile4);
WVEL_CE = load(infile5);
SSH_CE = load(infile7);
%===================


%% === data analysis ===
% convert lat and lon to x and y
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

% Lat_vec=U_CE.Lat_vec;Lon_vec=U_CE.Lon_vec;
% Lon=U_CE.Lon;Lat=U_CE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);


% zU_raw   = UVEL_CE.z; zU_raw(indxZ),
zW_raw   = WVEL_CE.z; zW_raw(indxZ) 

% u_ce_raw = squeeze(UVEL_CE.u(:,:,indxZ,:));
% v_ce_raw = squeeze(VVEL_CE.v(:,:,indxZ,:));

w_ce_raw = squeeze(WVEL_CE.w(:,:,indxZ,:));
eta_raw = SSH_CE.ssh;

taux_raw = TAUX_CE.taux;
tauy_raw = TAUY_CE.tauy;
tau_raw  = sqrt(taux_raw.^2+tauy_raw.^2);
tau_llav = squeeze(nanmean(nanmean(tau_raw,1),2)); % llav: lat and lon av

jultime_vec = TAUX_CE.jultime_vec;
jultime= datenum(jultime_vec);

clear TAUX_CE TAUY_CE UVEL_CE VVEL_CE WVEL_CE
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));
g = sw_g(Lat_vec(1),0); % g = sw_g(lat,z)

% compute geostrophic velocity and geostrophic vorticity
for i = 1 : length(jultime) % loop through time
    
    % compute geostrophic velocity
    [u_gs_raw(:,:,i),v_gs_raw(:,:,i)] = CalcGeostrophicVelFunc( ...
       x,y,eta_raw(:,:,i),f,g); 
    
    % compute geostrophic vorticity
    kesai_raw(:,:,i) = CalcRelVorticityFunc(x,y, ...
       u_gs_raw(:,:,i),v_gs_raw(:,:,i));
   
    Ro_raw(:,:,i) = kesai_raw(:,:,i)./f; % Rossby num
end

% compute linear and nonlinear Ekman Wvel
for i = 1 : length(jultime) % loop through time
    
    [w_nl_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
       taux_raw(:,:,i),tauy_raw(:,:,i),f,kesai_raw(:,:,i));
 
    [w_li_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
       taux_raw(:,:,i),tauy_raw(:,:,i),f,0);
 
end


figure
 pcolor(Lon,Lat,w_nl_raw(:,:,366));shading interp;colorbar;
 caxis([-5 5].*1e-5);polarmap  
