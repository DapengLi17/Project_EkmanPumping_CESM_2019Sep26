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

  date_str='2019Nov22';
  
infile1 = '../data_after_manipulation/CESM_TauUVWTw_dpthindxof5_KE_2019Nov22.mat';
% infile1 = '../data_after_manipulation/CESM_TauUVWTw_dpthindxof3_KE_2019Nov22.mat';

addpath(genpath('Func4EkmanProject/'))
rho_w = 1020;
%=============================


%% === load data ===
KE = load(infile1); % KE : Kuroshio Extension
%===================


%% === data analysis ===
% convert lat and lon to x and y
Lat_vec=KE.Lat_vec; Lon_vec=KE.Lon_vec;
Lon=KE.Lon; Lat=KE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);

u_ce_raw = KE.u;
v_ce_raw = KE.v;
w_ce_raw = KE.w;
taux = KE.taux;
tauy = KE.tauy;
tauA = sqrt(taux.^2+tauy.^2);
tauA_av = squeeze(nanmean(nanmean(tauA,1),2)); % whos tauA_av
Tw_raw = KE.temp;

jultime_vec= KE.jultime_vec;
jultime= datenum(jultime_vec);
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
for i = 1 : size(u_ce_raw,3) % loop through time
    kesai_ce(:,:,i) = CalcRelVorticityFunc(x,y,u_ce_raw(:,:,i),v_ce_raw(:,:,i));
end

% compute Rossby number
for i = 1 : size(kesai_ce,3)
    Ro_ce(:,:,i) = kesai_ce(:,:,i)./f;  
end

% compute nonlinear Ekman Wvel (Stern 1965)
for i = 1 : size(kesai_ce,3)
 [w_nl_raw(:,:,i)] = CalcEkmanWvelStern65Func(rho_w,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,kesai_ce(:,:,i));
end

% compute linear Ekman Wvel
for i = 1 : size(kesai_ce,3)
 [w_li_raw(:,:,i)] = CalcEkmanWvelStern65Func(rho_w,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,0);
end

w_df_raw = w_nl_raw - w_li_raw; % df: difference between high and low CESM model results

% find specific time (winter season)

  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{12}]; % Winter
  IndxSea{2} = [IndxMon{3};IndxMon{4};IndxMon{5}]; % Spring 
  IndxSea{3} = [IndxMon{6};IndxMon{7};IndxMon{8}]; % Summer
  IndxSea{4} = [IndxMon{9};IndxMon{10};IndxMon{11}]; % Fall
  
% whos w_nl  

w_nl = w_nl_raw(:,:,IndxSea{1});
w_ce = w_ce_raw(:,:,IndxSea{1});
w_df = w_df_raw(:,:,IndxSea{1});
w_li = w_li_raw(:,:,IndxSea{1});
Tw = Tw_raw(:,:,IndxSea{1});

JqEdy_nl = CalcVerticalEddyHeatFluxFunc(w_nl,Tw,rho_w);
JqEdy_nl_timeav = nanmean(JqEdy_nl,3);
JqEdy_nl_llav = squeeze(nanmean(nanmean(JqEdy_nl,1),2));

JqEdy_li = CalcVerticalEddyHeatFluxFunc(w_li,Tw,rho_w);
JqEdy_li_av = nanmean(JqEdy_li,3);
JqEdy_li_llav = squeeze(nanmean(nanmean(JqEdy_li,1),2));

JqEdy_df = CalcVerticalEddyHeatFluxFunc(w_df,Tw,rho_w);
JqEdy_df_av = nanmean(JqEdy_df,3);
JqEdy_df_llav = squeeze(nanmean(nanmean(JqEdy_df,1),2));

JqEdy_ce = CalcVerticalEddyHeatFluxFunc(w_ce,Tw,rho_w);
JqEdy_ce_av = nanmean(JqEdy_ce,3);
JqEdy_ce_llav = squeeze(nanmean(nanmean(JqEdy_ce,1),2));

tau_ampSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(tau_A,IndxSea{i},'nanmean');

figure;
 subplot(4,1,1)
 plotyy([1:90],JqEdy_nl_llav, ...
     [1:90],tauA_av(IndxSea{1}));
 subplot(4,1,2)
 plotyy([1:90],JqEdy_ce_llav, ...
     [1:90],tauA_av(IndxSea{1}));
 subplot(4,1,3)
 plot([1:90],JqEdy_ce_llav,'r');hold on;
 plot([1:90],JqEdy_nl_llav,'b'); 
 plot([1:90],JqEdy_li_llav,'k');
 plot([1:90],JqEdy_df_llav,'y');
 
 subplot(4,1,4)
 plot([1:90],JqEdy_nl_llav./JqEdy_ce_llav,'b');
 
 hold on;plot(tauA_av(IndxSea{1}))

nanmean(JqEdy_ce_av(:))
nanmean(JqEdy_nl_timeav(:))
nanmean(JqEdy_df_av(:))
nanmean(JqEdy_li_av(:))


indx = find(tauA_av>0.25); % length(indx) = 47
JqEdy_df_st = JqEdy_df(:,:,indx);
JqEdy_df_st_av = nanmean(JqEdy_df,3);

JqEdy_ce_st = JqEdy_ce(:,:,indx);
JqEdy_ce_st_av = nanmean(JqEdy_ce,3);


JqEdy_nl_st = JqEdy_nl(:,:,indx);
JqEdy_nl_st_av = nanmean(JqEdy_nl,3);

Jq_lim = [-100 100];

figure
  pcolor(Lon,Lat,JqEdy_ce_av);shading interp;
  caxis(Jq_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);

figure
  pcolor(Lon,Lat,JqEdy_df_av);shading interp;
  caxis(Jq_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);  
 
figure
  pcolor(Lon,Lat,JqEdy_nl_av);shading interp;
  caxis(Jq_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);   

figure
  pcolor(Lon,Lat,JqEdy_li_av);shading interp;
  caxis(Jq_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);    
  
% average during storm events
indx = find(tauA_av>0.2); % length(indx) = 47

w_li_st = w_li_raw(:,:,indx);
w_li_st_av = nanmean(w_li_st,3);
w_nl_st = w_nl_raw(:,:,indx);
w_nl_st_av = nanmean(w_nl_st,3);

w_df_st = w_df_raw(:,:,indx);
w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce_raw(:,:,indx);
w_ce_st_av = nanmean(w_ce_st,3);

u_ce_st = u_ce_raw(:,:,indx);
u_ce_st_av = nanmean(u_ce_st,3);

v_ce_st = v_ce_raw(:,:,indx);
v_ce_st_av = nanmean(v_ce_st,3);

Ro_ce_st = Ro_ce(:,:,indx);
Ro_ce_st_av = nanmean(Ro_ce_st,3);


plt_indx = 51;
figure
  pcolor(Lon,Lat,JqEdy_ce(:,:,plt_indx));shading interp;
  caxis(Jq_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);  
  
figure  
 pcolor(Lon,Lat,tauA(:,:,plt_indx));shading interp;
   caxis([0 0.25]);colorbar;colormap(ax1,parula);

% ================