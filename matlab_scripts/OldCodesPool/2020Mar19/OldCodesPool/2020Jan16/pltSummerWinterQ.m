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

  date_str='2019Dec05';
  
infile1 = '../data_after_manipulation/CESMtaux_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMtauy_KE_2019Dec01.mat';  
infile3 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
infile4 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';
infile5 = '../data_after_manipulation/CESMwvel_KE_2019Dec01.mat';
infile6 = '../data_after_manipulation/CESMtemp_KE_2019Dec01.mat';

pic1 = ['../pics/QekmanContours_KE_',date_str,'.png'];
pic2 = ['../pics/QekmanScatter_KE_',date_str,'.png'];
pic3 = ['../pics/QekmanPDF_KE_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))
addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))
rho_w = 1020;

% lat and lon range
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
Ro_std_limit = 0.15; 
tau_st_limt = 0.25; 
%=============================


%% === load data ===
TAUX_CE = load(infile1);
TAUY_CE = load(infile2);
UVEL_CE = load(infile3);
VVEL_CE = load(infile4);
WVEL_CE = load(infile5);
TEMP_CE = load(infile6);
%===================


%% === data analysis ===
% convert lat and lon to x and y
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

% Lat_vec=U_CE.Lat_vec;Lon_vec=U_CE.Lon_vec;
% Lon=U_CE.Lon;Lat=U_CE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);

zU_raw   = UVEL_CE.z;
zW_raw   = WVEL_CE.z; 
zT_raw   = TEMP_CE.z;
indxZ = 3;
zU_raw(indxZ),zW_raw(indxZ),zT_raw(indxZ),

u_ce_raw = squeeze(UVEL_CE.u(:,:,indxZ,:));
v_ce_raw = squeeze(VVEL_CE.v(:,:,indxZ,:));
w_ce_raw = squeeze(WVEL_CE.w(:,:,indxZ,:));
t_ce_raw = squeeze(TEMP_CE.temp(:,:,indxZ,:));

taux_raw = TAUX_CE.taux;
tauy_raw = TAUY_CE.tauy;
tau_raw  = sqrt(taux_raw.^2+tauy_raw.^2);
tau_llav = squeeze(nanmean(nanmean(tau_raw,1),2)); % llav: lat and lon av

jultime_vec = UVEL_CE.jultime_vec;
jultime= datenum(jultime_vec);

clear TAUX_CE TAUY_CE UVEL_CE VVEL_CE WVEL_CE TEMP_CE
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
for i = 1 : length(jultime) % loop through time
    kesai_raw(:,:,i) = CalcRelVorticityFunc(x,y, ...
       u_ce_raw(:,:,i),v_ce_raw(:,:,i));
   
    Ro_raw(:,:,i) = kesai_raw(:,:,i)./f; % Rossby num
end

% compute linear and nonlinear Ekman Wvel
for i = 1 : length(jultime) % loop through time
    
 [w_nl_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,kesai_raw(:,:,i));
 
 [w_li_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,0);
 
end

w_df_raw = w_nl_raw - w_li_raw; % df: difference between high and low CESM model

%%

% find specific time (winter season)
  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
      IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
  IndxSea{2} = [IndxMon{4};IndxMon{5};IndxMon{6}; ...
      IndxMon{7};IndxMon{8};IndxMon{9}]; % Winter season (Jing et al. 2019)

%   IndxSea{1} = [IndxMon{1};IndxMon{2}; ...
%       IndxMon{12}]; % Winter season (Jing et al. 2019)
  
% --- winter ---  
jultime_wt = jultime(IndxSea{1});

tau_llav_wt  = tau_llav(IndxSea{1});
tau_wt       = tau_raw(:,:,IndxSea{1});

kesai_wt = kesai_raw(:,:,IndxSea{1});
Ro_wt  = Ro_raw(:,:,IndxSea{1});

w_nl_wt = w_nl_raw(:,:,IndxSea{1});
w_li_wt = w_li_raw(:,:,IndxSea{1});
w_df_wt = w_df_raw(:,:,IndxSea{1});
w_ce_wt = w_ce_raw(:,:,IndxSea{1});
t_ce_wt = t_ce_raw(:,:,IndxSea{1});

[Q_li_wt,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_wt,t_ce_wt,rho_w);
% Q_nl_tmav_wt = nanmean(Q_nl_wt,3);
% Q_nl_llav_wt = squeeze(nanmean(nanmean(Q_nl_wt,1),2));


[Q_nl_wt,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_wt,t_ce_wt,rho_w);
% Q_nl_tmav_wt = nanmean(Q_nl_wt,3);
% Q_nl_llav_wt = squeeze(nanmean(nanmean(Q_nl_wt,1),2));

[Q_ce_wt,~,~]= CalcVerticalEddyHeatFluxFunc(w_ce_wt,t_ce_wt,rho_w);
% Q_ce_tmav_wt = nanmean(Q_ce_wt,3);
% Q_ce_llav_wt = squeeze(nanmean(nanmean(Q_ce_wt,1),2));
% ---------------

% --- summer ---  
jultime_sm = jultime(IndxSea{2});

tau_llav_sm  = tau_llav(IndxSea{2});
tau_sm       = tau_raw(:,:,IndxSea{2});

kesai_sm = kesai_raw(:,:,IndxSea{2});
Ro_sm  = Ro_raw(:,:,IndxSea{2});

w_nl_sm = w_nl_raw(:,:,IndxSea{2});
w_li_sm = w_li_raw(:,:,IndxSea{2});
w_df_sm = w_df_raw(:,:,IndxSea{2});
w_ce_sm = w_ce_raw(:,:,IndxSea{2});
t_ce_sm = t_ce_raw(:,:,IndxSea{2});

[Q_li_sm,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_sm,t_ce_sm,rho_w);
% Q_nl_tmav_sm = nanmean(Q_nl_sm,3);
% Q_nl_llav_sm = squeeze(nanmean(nanmean(Q_nl_sm,1),2));


[Q_nl_sm,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_sm,t_ce_sm,rho_w);
% Q_nl_tmav_sm = nanmean(Q_nl_sm,3);
% Q_nl_llav_sm = squeeze(nanmean(nanmean(Q_nl_sm,1),2));

[Q_ce_sm,~,~]= CalcVerticalEddyHeatFluxFunc(w_ce_sm,t_ce_sm,rho_w);
% Q_ce_tmav_sm = nanmean(Q_ce_sm,3);
% Q_ce_llav_sm = squeeze(nanmean(nanmean(Q_ce_sm,1),2));
% --------------------

%% === make pics ===
figure;

Q_lim_ce = [-400 400];
Q_lim_nl = [-200 200];

 subplot(2,3,1)
  pcolor(Lon,Lat,Q_ce_sm);shading interp;
    caxis(Q_lim_ce);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(2,3,2)
  pcolor(Lon,Lat,Q_ce_wt);shading interp;
    caxis(Q_lim_ce);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
    hold on;
  [c,h]=contour(Lon,Lat,Ro_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(2,3,3)
  pcolor(Lon,Lat,Q_ce_wt-Q_ce_sm);shading interp;
    caxis(Q_lim_ce);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(2,3,4)
  pcolor(Lon,Lat,Q_nl_sm);shading interp;
    caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(2,3,5)
  pcolor(Lon,Lat,Q_nl_wt);shading interp;
    caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(2,3,6)
  pcolor(Lon,Lat,Q_nl_wt-Q_nl_sm);shading interp;
    caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    
    
    Q_ce_df=Q_ce_wt-Q_ce_sm;% nanmean(Q_ce_df(:)) = 34.6360
    Q_nl_df=Q_nl_wt-Q_nl_sm;% nanmean(Q_nl_df(:)) = 1.0422
    
    Ro_std = nanstd(Ro_raw,0,3);
    
    indxRo = find(Ro_std>Ro_std_limit);length(indxRo)
    % nanmean(Q_ce_df(indxRo)) = 77.6909
    % nanmean(Q_nl_df(indxRo)) = 5.6599
    
% ==================
