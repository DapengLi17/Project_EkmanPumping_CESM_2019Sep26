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

% --- winter --- 
% find specific time (winter season)
  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
      IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
  
jultime_wt = jultime(IndxSea{1});

tau_llav_wt  = tau_llav(IndxSea{1});
tau_wt       = tau_raw(:,:,IndxSea{1});

kesai_wt = kesai_raw(:,:,IndxSea{1});
Ro_wt  = Ro_raw(:,:,IndxSea{1});
% Ro_st = Ro_wt(:,:,indx_st);
% Ro_st_av = nanmean(Ro_st,3);
Ro_wt_std = nanstd(Ro_wt,0,3);

w_nl_wt = w_nl_raw(:,:,IndxSea{1});
w_li_wt = w_li_raw(:,:,IndxSea{1});
w_df_wt = w_df_raw(:,:,IndxSea{1});
w_ce_wt = w_ce_raw(:,:,IndxSea{1});
t_ce_wt = t_ce_raw(:,:,IndxSea{1});  
  
% --- winter storm events ---
indx_st = find(tau_llav_wt>=tau_st_limt); % length(indx_st) = 353

tau_st = tau_wt(:,:,indx_st);
tau_st_av = nanmean(tau_st,3);

w_li_st = w_li_wt(:,:,indx_st);
w_li_st_av = nanmean(w_li_st,3);
w_nl_st = w_nl_wt(:,:,indx_st);
w_nl_st_av = nanmean(w_nl_st,3);

% w_df_st = w_df_wt(:,:,indx_st);
% w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce_wt(:,:,indx_st);
w_ce_st_av = nanmean(w_ce_st,3);

% kesai_st = kesai_wt(:,:,indx_st);
% kesai_std_st = nanstd(kesai_st,0,3);
t_ce_st = t_ce_wt(:,:,indx_st);

[Q_li_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_st,t_ce_st,rho_w);
% Q_nl_tmav_st = nanmean(Q_nl_st,3);
% Q_nl_llav_st = squeeze(nanmean(nanmean(Q_nl_st,1),2));

[Q_nl_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_st,t_ce_st,rho_w);
% Q_nl_tmav_st = nanmean(Q_nl_st,3);
% Q_nl_llav_st = squeeze(nanmean(nanmean(Q_nl_st,1),2));

[Q_ce_st,~,~]= CalcVerticalEddyHeatFluxFunc(w_ce_st,t_ce_st,rho_w);
% ----------------------

% --- non winter storm events ---
indx_ns = find(tau_llav_wt<tau_st_limt); % length(indx_ns) = 375

tau_ns = tau_wt(:,:,indx_ns);
tau_ns_av = nanmean(tau_ns,3);

w_li_ns = w_li_wt(:,:,indx_ns);
w_li_ns_av = nanmean(w_li_ns,3);
w_nl_ns = w_nl_wt(:,:,indx_ns);
w_nl_ns_av = nanmean(w_nl_ns,3);

% w_df_ns = w_df_wt(:,:,indx_ns);
% w_df_ns_av = nanmean(w_df_ns,3);
w_ce_ns = w_ce_wt(:,:,indx_ns);
w_ce_ns_av = nanmean(w_ce_ns,3);

t_ce_ns = t_ce_wt(:,:,indx_ns);

[Q_li_ns,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_ns,t_ce_ns,rho_w);
% Q_nl_tmav_ns = nanmean(Q_nl_ns,3);
% Q_nl_llav_ns = squeeze(nanmean(nanmean(Q_nl_ns,1),2));

[Q_nl_ns,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_ns,t_ce_ns,rho_w);
% Q_nl_tmav_ns = nanmean(Q_nl_ns,3);
% Q_nl_llav_ns = squeeze(nanmean(nanmean(Q_nl_ns,1),2));

[Q_ce_ns,~,~]= CalcVerticalEddyHeatFluxFunc(w_ce_ns,t_ce_ns,rho_w);
% ---------------------------------
    
% --- Ro > Ro std ---
indxRo = find(Ro_wt_std>Ro_std_limit); % length(indxRo) = 7152
% w_li_Ro = w_li_st_av(indxRo);
% w_nl_Ro = w_nl_st_av(indxRo); 
% w_ce_Ro = w_ce_st_av(indxRo); 
% Q_li_Ro = Q_li_st_av(indxRo);
% Q_nl_Ro = Q_nl_st_av(indxRo); 
% Q_ce_Ro = Q_ce_st_av(indxRo);  

    Q_ce_df=Q_ce_st-Q_ce_ns;% nanmean(Q_ce_df(:)) = 26.7935
    Q_nl_df=Q_nl_st-Q_nl_ns;% nanmean(Q_nl_df(:)) = 2.2183
    Q_li_df=Q_li_st-Q_li_ns;% nanmean(Q_li_df(:)) = 0.0241
    % nanmean(Q_ce_df(indxRo)) = 78.1108
    % nanmean(Q_nl_df(indxRo)) = 21.6455
    % nanmean(Q_li_df(indxRo)) = 6.9012

%%
figure;

Q_lim_ce = [-400 400];
Q_lim_nl = [-200 200];
Q_lim_li = [-100 100];

 subplot(3,3,1)
  pcolor(Lon,Lat,Q_ce_ns);shading interp;
    caxis(Q_lim_ce);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,2)
  pcolor(Lon,Lat,Q_ce_st);shading interp;
    caxis(Q_lim_ce);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,3)
  pcolor(Lon,Lat,Q_ce_st-Q_ce_ns);shading interp;
    caxis(Q_lim_ce);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,4)
  pcolor(Lon,Lat,Q_nl_ns);shading interp;
    caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,5)
  pcolor(Lon,Lat,Q_nl_st);shading interp;
    caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,6)
  pcolor(Lon,Lat,Q_nl_st-Q_nl_ns);shading interp;
    caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,7)
  pcolor(Lon,Lat,Q_li_ns);shading interp;
    caxis(Q_lim_li);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,8)
  pcolor(Lon,Lat,Q_li_st);shading interp;
    caxis(Q_lim_li);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
 subplot(3,3,9)
  pcolor(Lon,Lat,Q_li_st-Q_li_ns);shading interp;
    caxis(Q_lim_li);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);    
    
   

