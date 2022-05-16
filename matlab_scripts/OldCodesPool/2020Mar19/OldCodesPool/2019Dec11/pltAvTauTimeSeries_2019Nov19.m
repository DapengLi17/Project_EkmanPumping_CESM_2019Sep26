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
  
infile1 = '../data_after_manipulation/CESM_TauUVWTw_dpthindxof5_KE_2019Nov22.mat';

addpath(genpath('Func4EkmanProject/'))
rho0 = 1020;
%=============================


%% === load data ===
TAU = load(infile1); % KE : Kuroshio Extension
%===================


%% === data analysis ===
% convert lat and lon to x and y
Lat_vec=KE.Lat_vec;Lon_vec=KE.Lon_vec;
Lon=KE.Lon;Lat=KE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);

u_ce = KE.u;
v_ce = KE.v;
w_ce = KE.w;
taux = KE.taux;
tauy = KE.tauy;
tauA = sqrt(taux.^2+tauy.^2);
tauA_av = squeeze(nanmean(nanmean(tauA,1),2));
Tw = KE.temp;

jultime_vec= KE.jultime_vec;
jultime= datenum(jultime_vec);
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
for i = 1 : size(u_ce,3) % loop through time
    kesai_ce(:,:,i) = CalcRelVorticityFunc(x,y,u_ce(:,:,i),v_ce(:,:,i));
end

% compute Rossby number
for i = 1 : size(kesai_ce,3)
    Ro_ce(:,:,i) = kesai_ce(:,:,i)./f;  
end

% compute nonlinear Ekman Wvel (Stern 1965)
for i = 1 : size(kesai_ce,3)
 [w_nl(:,:,i)] = CalcEkmanWvelStern65Func(rho0,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,kesai_ce(:,:,i));
end

% compute linear Ekman Wvel
for i = 1 : size(kesai_ce,3)
 [w_li(:,:,i)] = CalcEkmanWvelStern65Func(rho0,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,0);
end

w_df = w_nl - w_li; % df: difference between high and low CESM model results

% average during storm events
indx = find(tauA_av>0.2); % length(indx) = 47

w_li_st = w_li(:,:,indx);
w_li_st_av = nanmean(w_li_st,3);
w_nl_st = w_nl(:,:,indx);
w_nl_st_av = nanmean(w_nl_st,3);

w_df_st = w_df(:,:,indx);
w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce(:,:,indx);
w_ce_st_av = nanmean(w_ce_st,3);

u_ce_st = u_ce(:,:,indx);
u_ce_st_av = nanmean(u_ce_st,3);

v_ce_st = v_ce(:,:,indx);
v_ce_st_av = nanmean(v_ce_st,3);

Ro_ce_st = Ro_ce(:,:,indx);
Ro_ce_st_av = nanmean(Ro_ce_st,3);
% ================


%% === make pics ===
f1=figure;
  plot(jultime,tauA_av,'b-');hold on;
  plot(jultime(indx),tauA_av(indx),'ro');grid on;
  datetick('x','dd/mmm','keepticks');
  
  W_lim = [-1e-5 1e-5];
  indxX = 1 : 5 : length(Lon_vec);
  indxY = 1 : 5 : length(Lat_vec);
  
  
f2=figure; % winter
  set(f2,'units','normalized','position',[0,0,1,1])
    
 subplot(3,2,1)
   quiver(Lon(indxY,indxX),Lat(indxY,indxX), ...
      u_ce_st_av(indxY,indxX),v_ce_st_av(indxY,indxX), ...
      'b','AutoScale','off');hold on;
   quiver(132,44,1,0,'b','AutoScale','off');
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
   text(134,44,'1 m/s');
   title(['CESM UV vel [m/s] dpth:',num2str(KE.dpth_u),'m']);
     
 subplot(3,2,2)
  pcolor(Lon,Lat,Ro_ce_st_av);shading interp;
  caxis([-0.3 0.3]);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['Rossby Num: kesai/f']);
  
 subplot(3,2,3)
  pcolor(Lon,Lat,w_li_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(linear) [m/s]']);
  
 subplot(3,2,4)
  pcolor(Lon,Lat,w_nl_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear) [m/s]']);
  
 subplot(3,2,5)
  pcolor(Lon,Lat,w_df_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear)-W(linear) [m/s]']);
  
 subplot(3,2,6)
  pcolor(Lon,Lat,w_ce_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);

 figure;
  pcolor(Lon,Lat,w_df_st_av./w_ce_st_av);shading interp;
  caxis([0 1]);
  
  w_ce_st_av_1d=w_ce_st_av(:);
  w_df_st_av_1d=w_df_st_av(:);
  m = find(isnan(w_ce_st_av_1d));
  n = find(isnan(w_df_st_av_1d));
  indxNaN = union(m,n);
  w_ce_st_av_1d(indxNaN)=[];
  w_df_st_av_1d(indxNaN)=[];  
  
  corr2(w_ce_st_av_1d,w_df_st_av_1d)

%% === make pics ===

indx_c = 18;

taux_st = taux(:,:,indx);
taux_st_av = taux_st(:,:,indx_c);

tauy_st = tauy(:,:,indx);
tauy_st_av = tauy_st(:,:,indx_c);

w_li_st = w_li(:,:,indx);
w_li_st_av = w_li_st(:,:,indx_c);
w_nl_st = w_nl(:,:,indx);
w_nl_st_av = w_nl_st(:,:,indx_c);

w_df_st = w_df(:,:,indx);
w_df_st_av = w_df_st(:,:,indx_c);

w_ce_st = w_ce(:,:,indx);
w_ce_st_av = w_ce_st(:,:,indx_c);

u_ce_st = u_ce(:,:,indx);
u_ce_st_av = u_ce_st(:,:,indx_c);

v_ce_st = v_ce(:,:,indx);
v_ce_st_av = v_ce_st(:,:,indx_c);

Ro_ce_st = Ro_ce(:,:,indx);
Ro_ce_st_av = Ro_ce_st(:,:,indx_c);


f1=figure;
quiver(Lon(indxY,indxX),Lat(indxY,indxX), ...
      taux_st_av(indxY,indxX),tauy_st_av(indxY,indxX), ...
      'b','AutoScale','off');hold on;
  
  W_lim = [-1e-5 1e-5];
  indxX = 1 : 5 : length(Lon_vec);
  indxY = 1 : 5 : length(Lat_vec);
  
  
f2=figure; % winter
  set(f2,'units','normalized','position',[0,0,1,1])
    
 subplot(3,2,1)
   quiver(Lon(indxY,indxX),Lat(indxY,indxX), ...
      u_ce_st_av(indxY,indxX),v_ce_st_av(indxY,indxX), ...
      'b','AutoScale','off');hold on;
   quiver(132,44,1,0,'b','AutoScale','off');
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
   text(134,44,'1 m/s');
   title(['CESM UV vel [m/s] dpth:',num2str(KE.dpth_u),'m']);
     
 subplot(3,2,2)
  pcolor(Lon,Lat,Ro_ce_st_av);shading interp;
  caxis([-0.3 0.3]);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['Rossby Num: kesai/f']);
  
 subplot(3,2,3)
  pcolor(Lon,Lat,w_li_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(linear) [m/s]']);
  
 subplot(3,2,4)
  pcolor(Lon,Lat,w_nl_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear) [m/s]']);
  
 subplot(3,2,5)
  pcolor(Lon,Lat,w_df_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear)-W(linear) [m/s]']);
  
 subplot(3,2,6)
  pcolor(Lon,Lat,w_ce_st_av);shading interp;
  caxis(W_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);

 figure;
  pcolor(Lon,Lat,w_df_st_av./w_ce_st_av);shading interp;
  caxis([0 1]);
  
  w_ce_st_av_1d=w_ce_st_av(:);
  w_df_st_av_1d=w_df_st_av(:);
  m = find(isnan(w_ce_st_av_1d));
  n = find(isnan(w_df_st_av_1d));
  indxNaN = union(m,n);
  w_ce_st_av_1d(indxNaN)=[];
  w_df_st_av_1d(indxNaN)=[];  
  
  corr2(w_ce_st_av_1d,w_df_st_av_1d)  