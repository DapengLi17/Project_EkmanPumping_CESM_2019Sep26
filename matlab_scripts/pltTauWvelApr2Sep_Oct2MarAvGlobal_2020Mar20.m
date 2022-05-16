%% === readme ===

% descrip: matlab scripts plot (1) wind stress amp and curl, 
% (2)linear, nonlinear Ekman W vel, and total W vel averaged over
% winter and summer season

% update history:
% v1.0 DL 2020Feb09

% extra notes:
% ===============


%% === set up environments ===
clear all;close all;clc;

  date_str='2020Mar20';

infile1 = '../data_after_manipulation/TauKesaiRoWvelTempGlobalAprSepAv_2020Mar17.mat';
infile2 = '../data_after_manipulation/TauKesaiRoWvelTempGlobalOctMarAv_2020Mar17.mat';

pic1 = ['../pics/TauAmpCurlApr2Sep_Oct2MarAvGlobal_',date_str,'.png'];
pic2 = ['../pics/WvelEkmanCESMApr2Sep_Oct2MarAvGlobal_',date_str,'.png'];

  addpath(genpath('Func4EkmanProject/'))
  
% Constants4CESM_Global

 lon_limits = [0 360];
 lon_ticks = [0:60:360];
 lat_limits = [-70 70];
 lat_ticks = [-60:30:60];
%=============================


%% === load data ===
Apr2Sep = load(infile1);
Oct2Mar = load(infile2);
%===================


%% === data analysis ===
lat_1d_raw = Apr2Sep.lat_1d;
lon_1d_raw = Apr2Sep.lon_1d;

indxLat = find(lat_1d_raw >= lat_limits(1) & lat_1d_raw <= lat_limits(2));
lat_1d = lat_1d_raw(indxLat);
indxLon = find(lon_1d_raw >= lon_limits(1) & lon_1d_raw <= lon_limits(2));
lon_1d = lon_1d_raw(indxLon);

tau_amp_Apr2Sep = Apr2Sep.tau_amp_TimeAv(indxLat,indxLon);
tau_amp_Oct2Mar = Oct2Mar.tau_amp_TimeAv(indxLat,indxLon);
tau_curl_Apr2Sep = Apr2Sep.tau_curl_TimeAv(indxLat,indxLon);
tau_curl_Oct2Mar = Oct2Mar.tau_curl_TimeAv(indxLat,indxLon);

W_ce_Oct2Mar = Oct2Mar.w_ce_TimeAv(indxLat,indxLon);
W_ce_Apr2Sep = Apr2Sep.w_ce_TimeAv(indxLat,indxLon);
W_te_Oct2Mar = Oct2Mar.w_te_TimeAv(indxLat,indxLon);
W_te_Apr2Sep = Apr2Sep.w_te_TimeAv(indxLat,indxLon);
W_le_Oct2Mar = Oct2Mar.w_le_TimeAv(indxLat,indxLon);
W_le_Apr2Sep = Apr2Sep.w_le_TimeAv(indxLat,indxLon);
% W_ne_Oct2Mar = Oct2Mar.w_ne_TimeAv(indxLat,indxLon);
% W_ne_Apr2Sep = Apr2Sep.w_ne_TimeAv(indxLat,indxLon);
W_ne_Oct2Mar = W_te_Oct2Mar - W_le_Oct2Mar;
W_ne_Apr2Sep = W_te_Apr2Sep - W_le_Apr2Sep;

clear Apr2Sep Oct2Mar

% --- NaN the Equatorial Regions ---
indx4NaN = find(abs(lat_1d)<=10);

% size W_te_Oct2Mar: 1400x3600, lon_1d: 1x3600, lat_1d: 1400x1
W_te_Apr2Sep(indx4NaN,:)=NaN; W_te_Oct2Mar(indx4NaN,:)=NaN;
W_le_Apr2Sep(indx4NaN,:)=NaN; W_le_Oct2Mar(indx4NaN,:)=NaN;
W_ne_Apr2Sep(indx4NaN,:)=NaN; W_ne_Oct2Mar(indx4NaN,:)=NaN;
% ======================


%% === make pics ===
f1=figure('Renderer','painters','visible','off'); % wind stress amp

  pic1_size = [8.5 4]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
  cbar_limits_TauAmp  = [0 0.4];
  cbar_ticks_TauAmp  = [0:0.1:0.4];
  cbar_limits_TauCurl = [-3 3];
  cbar_ticks_TauCurl = [-3:1:3];
  x4text = -2.9;
  y4text = -1.3;
  ylabel_position=[-0.1, 0.5, 0];
  font_size = 8;
  
% ~~~ generate subplot position ~~~  
  row_num=2;col_num=2;
  margin_left=0.08;margin_right=0.08;
  margin_top=0.05;margin_botm=0.08;
  pics_dist_x=0.05; pics_dist_y=0.05;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

m_proj('miller','lon',lon_limits,'lat',lat_limits);

 subplot('Position',sbplt_posit(1,:))
   m_pcolor(lon_1d,lat_1d,tau_amp_Apr2Sep);shading interp;polarmap;
    caxis(cbar_limits_TauAmp);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Apr-Sep Average');
    ylabel('\tau','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,tau_amp_Oct2Mar);shading interp;
    caxis(cbar_limits_TauAmp);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Oct-Mar Average');
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_TauAmp,'YTick',cbar_ticks_TauAmp);
   title(h1,'[N/m^2]','FontSize',font_size-2);    
   
   
 subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,tau_curl_Apr2Sep.*10^7);shading interp;
    caxis(cbar_limits_TauCurl);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Curl(\tau)','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,tau_curl_Oct2Mar.*10^7);shading interp;
    caxis(cbar_limits_TauCurl);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
    h1=colorbar('FontSize',font_size-2);
    set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_TauCurl,'YTick',cbar_ticks_TauCurl);
    title(h1,'[x10^{-7} N/m^3]','FontSize',font_size-2);  

   set(gcf,'color','w');
% ========================


%% === make pics ===
f2=figure('Renderer','painters','visible','off'); % Wvel

 pic2_size = [8.5 6.5]; 
 set(f2,'units','inches','position',[5,5,pic2_size])
 cbar_limits_Wvel = [-0.5 0.5];
 cbar_ticks_Wvel = [-0.5:0.25:0.5];
 x4text = -2.9;
 y4text = -1.3;
 ylabel_position=[-0.1, 0.5, 0];
 font_size = 8;
  
 % ~~~ generate subplot position ~~~  
  row_num=4;col_num=2;
  margin_left=0.08;margin_right=0.08;
  margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.05; pics_dist_y=0.04;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

m_proj('miller','lon',lon_limits,'lat',lat_limits);

  subplot('Position',sbplt_posit(1,:))
   m_pcolor(lon_1d,lat_1d,W_ce_Apr2Sep.*86400);shading interp;polarmap;
    caxis(cbar_limits_Wvel.*2);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Apr-Sep Average');
    ylabel('W_{total}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,W_ce_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wvel.*2);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Oct-Mar Average');
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Wvel.*2,'YTick',cbar_ticks_Wvel.*2);
   title(h1,'[m/day]','FontSize',font_size-2);  
   
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,W_te_Apr2Sep.*86400);shading interp;
    caxis(cbar_limits_Wvel);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{TE}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,W_te_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wvel);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_Wvel,'YTick',cbar_ticks_Wvel);
   title(h1,'[m/day]','FontSize',font_size-2);  
   
   
  subplot('Position',sbplt_posit(5,:))
   m_pcolor(lon_1d,lat_1d,W_le_Apr2Sep.*86400);shading interp;
    caxis(cbar_limits_Wvel);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{LE}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'e','fontsize',font_size,'FontWeight','bold');
    
  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,W_le_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wvel);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'f','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',cbar_limits_Wvel,'YTick',cbar_ticks_Wvel);
   title(h1,'[m/day]','FontSize',font_size-2);  
   

  subplot('Position',sbplt_posit(7,:))   
   m_pcolor(lon_1d,lat_1d,W_ne_Apr2Sep.*86400);shading interp;
    caxis(cbar_limits_Wvel);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{NE}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'g','fontsize',font_size,'FontWeight','bold');
  
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,W_ne_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wvel);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'h','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(8,2) .01 sbplt_posit(8,4)-0.02], ...
       'ylim',cbar_limits_Wvel,'YTick',cbar_ticks_Wvel);
   title(h1,'[m/day]','FontSize',font_size-2);  

   set(gcf,'color','w');
% ====================

  
%% === output data ===   
export_fig(f1,pic1,'-r200','-nocrop')
export_fig(f2,pic2,'-r200','-nocrop')
% ==================== 
