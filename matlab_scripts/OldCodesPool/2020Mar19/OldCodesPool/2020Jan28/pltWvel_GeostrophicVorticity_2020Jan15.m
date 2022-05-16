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

w_df_raw = w_nl_raw - w_li_raw; % df: difference between high and low CESM model results

% find specific time (winter season)
  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% % Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
      IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
%   IndxSea{1} = [IndxMon{1};IndxMon{2}; ...
%       IndxMon{12}]; % Winter season
  
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
% ---------------

% composite av during winter storm events
indx_st = find(tau_llav_wt>=tau_st_limt); % length(indx_st) = 146
indx_ns = find(tau_llav_wt<tau_st_limt); % length(indx_ns) = 214

tau_st = tau_wt(:,:,indx_st);
tau_st_av = nanmean(tau_st,3);

Ro_st = Ro_wt(:,:,indx_st);
Ro_st_av = nanmean(Ro_st,3);
Ro_st_std = nanstd(Ro_st,0,3);

w_li_st = w_li_wt(:,:,indx_st);
w_li_st_av = nanmean(w_li_st,3);
w_nl_st = w_nl_wt(:,:,indx_st);
w_nl_st_av = nanmean(w_nl_st,3);

w_df_st = w_df_wt(:,:,indx_st);
w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce_wt(:,:,indx_st);
w_ce_st_av = nanmean(w_ce_st,3);

kesai_st = kesai_wt(:,:,indx_st);
kesai_std_st = nanstd(kesai_st,0,3); 
% ==================


%% === make pics ===
f1 = figure; 

  pic1_size=[12 4]; 
% pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  font_size=12;
  x_lim = [jultime(1) jultime(end)];
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
% 'Position',[5,5,pic1_size] unit is inch, 5,5, is the left-bottom corner
% of the fig on screen(no need to chage), pic_size specify the size of the window 
% as the same size of the output fig. 

% ~~~ generate subplot position ~~~  
  row_num=1;col_num=1;margin_left=0.06;
  margin_right=0.05;margin_top=0.1;margin_botm=0.1;
  pics_dist_x=0.0; pics_dist_y=0.0;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

subplot('Position',sbplt_posit(1,:));
 plot(jultime,tau_llav,'b','linewidth',1.5);hold on;
 plot(jultime([1,end]),[tau_st_limt tau_st_limt],'r--','linewidth',1.5);grid on;
 plot(jultime_wt(indx_st),tau_llav_wt(indx_st),'ro', ...
     'markersize',6,'MarkerFaceColor','r');
 plot(jultime_wt(indx_ns),tau_llav_wt(indx_ns),'go', ...
     'markersize',6,'MarkerFaceColor','g');
  set(gca,'xlim',x_lim,'fontsize',font_size);
  title(['spatial mean wind stress, ',num2str(length(indx_st)), ...
      ' winter storms totally']);ylabel('[N/m^2]','fontsize',font_size);
  datetick('x','mmm/yy');

%%
f3 = figure; 

  pic3_size=[4 4]; 
% pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  font_size=12;
  set(f3,'Units','inches','Position',[5,5,pic3_size]);
% 'Position',[5,5,pic1_size] unit is inch, 5,5, is the left-bottom corner
% of the fig on screen(no need to chage), pic_size specify the size of the window 
% as the same size of the output fig. 
  load coastlines

% ~~~ generate subplot position ~~~  
  row_num=1;col_num=1;margin_left=0.15;
  margin_right=0.08;margin_top=0.1;margin_botm=0.15;
  pics_dist_x=0.0; pics_dist_y=0.0;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

subplot('Position',sbplt_posit(1,:));
 plot(coastlon, coastlat,'k','linewidth',1.5);grid on;
 set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'ytick',[30:5:45])
 xlabel('Lon','fontsize',font_size);ylabel('Lat','fontsize',font_size);
 
 
%% 
f2=figure; 

  pic2_size=[10 10]; 
% pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  font_size=12;
  set(f2,'Units','inches','Position',[5,5,pic2_size]);
% 'Position',[5,5,pic1_size] unit is inch, 5,5, is the left-bottom corner
% of the fig on screen(no need to chage), pic_size specify the size of the window 
% as the same size of the output fig.
  W_lim = [-2e-5 2e-5];

% ~~~ generate subplot position ~~~  
  row_num=3;col_num=2;margin_left=0.06;
  margin_right=0.03;margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.03; pics_dist_y=0.05;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 
ax1=subplot('Position',sbplt_posit(1,:));
  pcolor(Lon,Lat,tau_st_av);shading interp;
  caxis([0 0.5]);colormap(ax1,polarmap(64,1));hc=colorbar;title(hc,'N/m^2');
  set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5],'fontsize',font_size);
  ylabel('Lat','fontsize',font_size);
  title([num2str(length(indx_st)),' storms composite-av \tau']);

ax2=subplot('Position',sbplt_posit(2,:));
  pcolor(Lon,Lat,Ro_st_av);shading interp;
  caxis([-0.3 0.3]);colormap(ax2,polarmap(64,1));
    hcb=colorbar;set(hcb,'ylim',[-0.3 0.3],'YTick',[-0.3:0.1:0.3])
    hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
  set(h,'LineWidth',1.5);
  set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'xticklabel',[], ...
      'yticklabel',[],'fontsize',font_size);
  title(['Ro']);

ax3=subplot('Position',sbplt_posit(3,:));
  pcolor(Lon,Lat,w_li_st_av);shading interp;
  caxis(W_lim);colormap(ax3,polarmap(64,1));hc=colorbar;title(hc,'m/s');
%     hold on;
%   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
%   set(h,'LineWidth',1.5);
  set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5],'fontsize',font_size);
    ylabel('Lat','fontsize',font_size);title(['W_f']);  
  
ax4=subplot('Position',sbplt_posit(4,:));
  pcolor(Lon,Lat,w_nl_st_av);shading interp;
%     hold on;
%   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
%   set(h,'LineWidth',1.5);
  caxis(W_lim);colormap(ax4,polarmap(64,1));hc=colorbar;title(hc,'m/s');
  set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5],'yticklabel',[],'fontsize',font_size);
  title(['W_{\xi}']);
  
ax5=subplot('Position',sbplt_posit(5,:));
  pcolor(Lon,Lat,w_df_st_av);shading interp;
    hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
  set(h,'LineWidth',1.5);
  caxis(W_lim);colormap(ax5,polarmap(64,1));hc=colorbar;title(hc,'m/s');
  set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'fontsize',font_size);
    xlabel('Lon','fontsize',font_size);ylabel('Lat','fontsize',font_size);
    title(['W_{\xi}-W_f']);
  
ax6=subplot('Position',sbplt_posit(6,:));
  pcolor(Lon,Lat,w_ce_st_av);shading interp;
    hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
   set(h,'LineWidth',1.5);
  caxis(W_lim);colormap(ax6,polarmap(64,1));hc=colorbar;title(hc,'m/s'); 
  set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'yticklabel',[],'fontsize',font_size);
  xlabel('Lon','fontsize',font_size);
  title(['W_{cesm} depth:',num2str(zW_raw(indxZ)),'m']); 
% ====================

  
%% === output data ===
print(f1,'-dpng','-r500',pic1) 
print(f2,'-dpng','-r500',pic2) 
print(f3,'-dpng','-r500',pic3) 
% ====================