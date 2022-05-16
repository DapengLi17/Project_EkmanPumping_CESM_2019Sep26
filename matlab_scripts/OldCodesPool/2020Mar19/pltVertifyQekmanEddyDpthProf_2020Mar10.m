%% === readme ===

% descrip: matlab scripts plot Qeddy at 50m from (1) one singer layer calculation, 
% (2) depth prof at 50m for comparison. 

% update history:
% v1.0 DL 2020Mar10

% extra notes:
% ===============


%% === set up environments ===
clear all;close all;clc;

  date_str = '2020Mar10';

infile1 = '../data_after_manipulation/QekmanEddyGlobalDpthProfOct-MarAv_2020Feb11.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalBorealWinterAv_2020Feb11.mat';
infile3 = '../data_after_manipulation/QekmanEddyGlobalDpthProfApr-SepAv_2020Feb11.mat';
infile4 = '../data_after_manipulation/QekmanEddyGlobalBorealSummerAv_2020Feb11.mat';

  pic1 = ['../pics/VerifyQeddyDpthProf7_',date_str,'.png'];
  pic2 = ['../pics/VerifyQekmanDpthProf7_',date_str,'.png'];
  
addpath(genpath('Func4EkmanProject/'))

  lon_limits = [0 360];
  lat_limits = [-70 70];
% =============================


%% === load data ===
DpthProf_wt = load(infile1);
Dpth50m_wt = load(infile2);
DpthProf_su = load(infile3);
Dpth50m_su = load(infile4);
% ==================


%% === data analysis ===
lon_1d = Dpth50m_wt.lon_1d;
lat_1d = Dpth50m_wt.lat_1d;
% DpthProf_wt.z_wvel_ce = 0 10 20 30 40 50 70 90 110
% DpthProf_wt.z_Tw = 5    15    25    35    45    55    75    95   115
Q_ce_prof50m_wt = DpthProf_wt.Q_ce(:,:,7)';
Q_ce_prof50m_su = DpthProf_su.Q_ce(:,:,7)';
Q_ce_dpth50m_wt = Dpth50m_wt.Q_ce;
Q_ce_dpth50m_su = Dpth50m_su.Q_ce;

Q_nl_prof50m_wt = DpthProf_wt.Q_nl(:,:,7)';
Q_nl_prof50m_su = DpthProf_su.Q_nl(:,:,7)';
Q_nl_dpth50m_wt = Dpth50m_wt.Q_nl;
Q_nl_dpth50m_su = Dpth50m_su.Q_nl;
% ======================


%% === make pics ===  
f1=figure('Renderer','painters','visible','off'); % Q

  pic1_size = [10 15]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
%   set(f1,'units','normalized','position',[0,0,1,1])
  font_size = 10;
  cbar_limits_Qce = [-100 100];
  cbar_ticks_Qce = [-100:50:100];
  cbar_limits_Qek = [-20 20]; 
  cbar_ticks_Qek = [-20:10:20];
  x4text = -2.9;
  y4text = -1.3;
  ylabel_position=[-0.08, 0.5, 0];
 
 % ~~~ generate subplot position ~~~  
  row_num=4; col_num=2;
  margin_left=0.08;margin_right=0.08;
  margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.05; pics_dist_y=0.04;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

m_proj('miller','lon',lon_limits,'lat',lat_limits);
   
  subplot('Position',sbplt_posit(1,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_dpth50m_su);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
    title('summer');
    ylabel('Q_{eddy}^{dpth:50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);

  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_dpth50m_wt);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
    title('winter');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Qce,'YTick',cbar_ticks_Qce);
   title(h1,'[W/m^2]','FontSize',font_size); 
   
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_prof50m_su);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
    ylabel('Q_{eddy}^{prof:50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
  
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_prof50m_wt);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_Qce,'YTick',cbar_ticks_Qce);
   title(h1,'[W/m^2]','FontSize',font_size);
   
  subplot('Position',sbplt_posit(5,:))
   m_pcolor(lon_1d,lat_1d,Q_nl_dpth50m_su);shading interp;
    polarmap;caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
    ylabel('Q_{nl-ek}^{dpth:50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);

  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,Q_nl_dpth50m_wt);shading interp;
    polarmap;caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size); 
   
  subplot('Position',sbplt_posit(7,:))
   m_pcolor(lon_1d,lat_1d,Q_nl_prof50m_su);shading interp;
    polarmap;caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
    ylabel('Q_{nl-ek}^{prof:50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
  
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,Q_nl_prof50m_wt);shading interp;
    polarmap;caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size);axis normal;
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(8,2) .01 sbplt_posit(8,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size); 
% =========================


%% === make pics ===
export_fig(f1,pic1,'-r200','-nocrop')    
% ==================