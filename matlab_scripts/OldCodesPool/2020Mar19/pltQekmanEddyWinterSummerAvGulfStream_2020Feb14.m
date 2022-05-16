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

   date_str='2020Feb14';
   
infile1 = '../data_after_manipulation/QekmanEddyGlobalBorealSummerAv_2020Feb11.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalBorealWinterAv_2020Feb11.mat';

outfile1 = ['../data_after_manipulation/QekmanEddyWinterSummerAvGulfStream_',date_str,'.txt'];
pic1 = ['../pics/QekmanEddyWinterSummerAvGulfStream_',date_str,'.png'];

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

   lat_limits = [25 50];  
   lon_limits = [275 325];
   lat_ticks = [30:10:50];
   lon_ticks = [280:10:320];
%=============================


%% === load data ===
Su = load(infile1);
Wt = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = Su.lat_1d;
lon_1d_raw = Su.lon_1d;

indxLat = find(lat_1d_raw >= lat_limits(1) & lat_1d_raw <= lat_limits(2));
lat_1d = lat_1d_raw(indxLat);
indxLon = find(lon_1d_raw >= lon_limits(1) & lon_1d_raw <= lon_limits(2));
lon_1d = lon_1d_raw(indxLon);

Q_li_wt = Wt.Q_li(indxLat,indxLon); 
Q_li_su = Su.Q_li(indxLat,indxLon); 
Q_nl_wt = Wt.Q_nl(indxLat,indxLon);
Q_nl_su = Su.Q_nl(indxLat,indxLon);
Q_ce_wt = Wt.Q_ce(indxLat,indxLon);
Q_ce_su = Su.Q_ce(indxLat,indxLon);

clear Wt Su
% ======================


%% === make pics ===
f1=figure('Renderer','painters'); % Q

  pic1_size = [8.5 6.5]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
%   set(f1,'units','normalized','position',[0,0,1,1])
  font_size = 8;
  cbar_limits_Qce = [-120 120];
  cbar_ticks_Qce = [-120:60:120];
  cbar_limits_Qek = [-40 40]; 
  cbar_ticks_Qek = [-40:20:40];
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
   m_pcolor(lon_1d,lat_1d,Q_ce_su);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
% --- save codes below as reference ---
%     set(gca,'color',[0.7,0.7,0.7], ... % set ax background color as gray for NaN (land)
%       'xlim',lon_limits,'xtick',lon_ticks,'XTickLabel',lon_tickslabel, ...
%       'ylim',lat_limits,'ytick',lat_ticks,'YTickLabel',lat_tickslabel, ...
%       'Fontsize',font_size);
% -------------------------------------
    title('Summer');
    ylabel('Q_{eddy}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_wt);shading interp;
    caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Winter');
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Qce,'YTick',cbar_ticks_Qce);
   title(h1,'[W/m^2]','FontSize',font_size-2);  

      
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,Q_nl_su);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{\xi}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,Q_nl_wt);shading interp;
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
   m_pcolor(lon_1d,lat_1d,Q_li_su);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{f}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'e','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,Q_li_wt);shading interp;
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
   m_pcolor(lon_1d,lat_1d,(Q_nl_su-Q_li_su));shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{\xi} - Q_{f}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'g','fontsize',font_size,'FontWeight','bold');
      
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,(Q_nl_wt-Q_li_wt));shading interp;
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
% ====================

  
%% === output data ===  
export_fig(f1,pic1,'-r200','-nocrop')

% write statistics to data file
fid = fopen(outfile1,'w');

 fprintf(fid,['readme: generated by DL on', date_str, ...
     ' via code pltQekmanEddyWinterSummerAvGulfStream_2020Feb14.m \n']);
 fprintf(fid,' \n');

 fprintf(fid,['                       Summer(Apr-Sep)  Winter(Oct-Mar) \n']);
 fprintf(fid,['                        mean+-std        mean+-std \n']);
 
%
 fprintf(fid,'Q_LIEkman(50m) [W/m2]   %.2f +- %.2f     %.2f +- %.2f\n', ...
   [nanmean(Q_li_su(:)),nanstd(Q_li_su(:)),nanmean(Q_li_wt(:)),nanstd(Q_li_wt(:))]);
 fprintf(fid,' \n');

%
 fprintf(fid,'Q_NLEkman(50m) [W/m2]   %.2f +- %.2f     %.2f +- %.2f\n', ...
   [nanmean(Q_nl_su(:)),nanstd(Q_nl_su(:)),nanmean(Q_nl_wt(:)),nanstd(Q_nl_wt(:))]);
 fprintf(fid,' \n'); 
 
%
 fprintf(fid,'Q_eddy(50m)    [W/m2]   %.2f +- %.2f   %.2f +- %.2f\n', ...
   [nanmean(Q_ce_su(:)),nanstd(Q_ce_su(:)),nanmean(Q_ce_wt(:)),nanstd(Q_ce_wt(:))]);
 fprintf(fid,' \n');
 
fclose(fid);  
% ====================
