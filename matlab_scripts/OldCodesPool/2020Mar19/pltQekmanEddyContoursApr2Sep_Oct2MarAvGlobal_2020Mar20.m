%% ===readme===

% descrip: matlab scripts plot W vel during winter storms from 
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

   date_str='2020Mar20';
   
infile1 = '../data_after_manipulation/QekmanEddyGlobalAprSepAv_2020Mar18.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalOctMarAv_2020Mar18.mat';

pic1 = ['../pics/QekmanEddyApr2Sep_Oct2MarAvGlobal_',date_str,'.png'];

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

 lon_limits = [0 360];
 lon_ticks = [0:60:360];
 lat_limits = [-70 70];
 lat_ticks = [-60:30:60];
%=============================


%% === load data ===
Apr2Sep = load(infile1);
Oct2Mar = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = Apr2Sep.lat_1d;
lon_1d_raw = Apr2Sep.lon_1d;

indxLat = find(lat_1d_raw >= lat_limits(1) & lat_1d_raw <= lat_limits(2));
lat_1d = lat_1d_raw(indxLat);
indxLon = find(lon_1d_raw >= lon_limits(1) & lon_1d_raw <= lon_limits(2));
lon_1d = lon_1d_raw(indxLon);

Q_te_Oct2Mar = Oct2Mar.Q_te(indxLat,indxLon);
Q_te_Apr2Sep = Apr2Sep.Q_te(indxLat,indxLon);
Q_le_Oct2Mar = Oct2Mar.Q_le(indxLat,indxLon);
Q_le_Apr2Sep = Apr2Sep.Q_le(indxLat,indxLon);
Q_ne_Oct2Mar = Oct2Mar.Q_ne(indxLat,indxLon);
Q_ne_Apr2Sep = Apr2Sep.Q_ne(indxLat,indxLon);
Q_ce_Oct2Mar = Oct2Mar.Q_ce(indxLat,indxLon);
Q_ce_Apr2Sep = Apr2Sep.Q_ce(indxLat,indxLon);

clear Apr2Sep Oct2Mar

% --- NaN the Equatorial Regions ---
indx4NaN = find(abs(lat_1d)<=10);

% size W_te_Oct2Mar: 1400x3600, lon_1d: 1x3600, lat_1d: 1400x1
Q_te_Apr2Sep(indx4NaN,:)=NaN; Q_te_Oct2Mar(indx4NaN,:)=NaN;
Q_le_Apr2Sep(indx4NaN,:)=NaN; Q_le_Oct2Mar(indx4NaN,:)=NaN;
Q_ne_Apr2Sep(indx4NaN,:)=NaN; Q_ne_Oct2Mar(indx4NaN,:)=NaN;
% ======================


%% === make pics ===
f1=figure('Renderer','painters','visible','off');

  pic1_size = [8.5 6.5]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
%   set(f1,'units','normalized','position',[0,0,1,1])
  font_size = 8;
  cbar_limits_Qce = [-100 100];
  cbar_ticks_Qce = [-100:50:100];
  cbar_limits_Qek = [-20 20]; 
  cbar_ticks_Qek = [-20:10:20];
  x4text = -2.9;
  y4text = -1.3;
  ylabel_position=[-0.1, 0.5, 0];
 
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
   m_pcolor(lon_1d,lat_1d,Q_ce_Apr2Sep);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Apr-Sep');
    ylabel('Q_{total}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_Oct2Mar);shading interp;
    caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Oct-Mar');
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Qce,'YTick',cbar_ticks_Qce);
   title(h1,'[W/m^2]','FontSize',font_size-2);  

      
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,Q_te_Apr2Sep);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{TE}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,Q_te_Oct2Mar);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size-2);  


  subplot('Position',sbplt_posit(5,:))
   m_pcolor(lon_1d,lat_1d,Q_le_Apr2Sep);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{LE}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'e','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,Q_le_Oct2Mar);shading interp;
    caxis(cbar_limits_Qek); 
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'f','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size-2);     
   
  subplot('Position',sbplt_posit(7,:))
   m_pcolor(lon_1d,lat_1d,Q_ne_Apr2Sep);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{NE}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'g','fontsize',font_size,'FontWeight','bold');
      
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,Q_ne_Oct2Mar);shading interp;
    caxis(cbar_limits_Qek); 
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'h','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(8,2) .01 sbplt_posit(8,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size-2);   

   set(gcf,'color','w');
% ====================

  
%% === output data ===  
export_fig(f1,pic1,'-r200','-nocrop')
% ====================
