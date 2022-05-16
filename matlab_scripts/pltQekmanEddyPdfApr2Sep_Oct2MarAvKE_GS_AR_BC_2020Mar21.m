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

   date_str='2020Mar21';
   
infile1 = '../data_after_manipulation/QekmanEddyGlobalAprSepAv_2020Mar18.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalOctMarAv_2020Mar18.mat';

pic1 = ['../pics/QekmanEddyPdfApr2Sep_Oct2MarAvKE_GS_AR_BC_',date_str,'.png'];

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

QbinsCenter = [-10:10:40]; % for PDF

% --- Kuroshio extension region ---
   lat_limits_KE = [30 45]; 
   lon_limits_KE = [140 170];
% ----------------------------------

% --- Gulf stream region ---
   lat_limits_GS = [35 42];  
   lon_limits_GS = [285 305];
% ---------------------------

% --- Agulhas return currents ---
   lat_limits_AR = [-45 -35]; 
   lon_limits_AR = [10 60];
% -------------------------------

% --- Brazil Current region ---
   lat_limits_BC = [-50 -35];  
   lon_limits_BC = [302 315];
% -----------------------------
%=============================


%% === load data ===
Apr2Sep = load(infile1);
Oct2Mar = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = Apr2Sep.lat_1d;
lon_1d_raw = Apr2Sep.lon_1d;

z = Apr2Sep.z_wvel_ce; 

indxLat_KE = find(lat_1d_raw >= lat_limits_KE(1) & lat_1d_raw <= lat_limits_KE(2));
indxLon_KE = find(lon_1d_raw >= lon_limits_KE(1) & lon_1d_raw <= lon_limits_KE(2));

indxLat_GS = find(lat_1d_raw >= lat_limits_GS(1) & lat_1d_raw <= lat_limits_GS(2));
indxLon_GS = find(lon_1d_raw >= lon_limits_GS(1) & lon_1d_raw <= lon_limits_GS(2));

indxLat_AR = find(lat_1d_raw >= lat_limits_AR(1) & lat_1d_raw <= lat_limits_AR(2));
indxLon_AR = find(lon_1d_raw >= lon_limits_AR(1) & lon_1d_raw <= lon_limits_AR(2));

indxLat_BC = find(lat_1d_raw >= lat_limits_BC(1) & lat_1d_raw <= lat_limits_BC(2));
indxLon_BC = find(lon_1d_raw >= lon_limits_BC(1) & lon_1d_raw <= lon_limits_BC(2));

Q_ce_Oct2Mar_KE = Oct2Mar.Q_ce(indxLat_KE,indxLon_KE);
Q_ce_Apr2Sep_KE = Apr2Sep.Q_ce(indxLat_KE,indxLon_KE);
Q_te_Oct2Mar_KE = Oct2Mar.Q_te(indxLat_KE,indxLon_KE);
Q_te_Apr2Sep_KE = Apr2Sep.Q_te(indxLat_KE,indxLon_KE);
% Q_le_Oct2Mar_KE = Oct2Mar.Q_le(indxLat_KE,indxLon_KE); 
% Q_le_Apr2Sep_KE = Apr2Sep.Q_le(indxLat_KE,indxLon_KE); 
% Q_ne_Oct2Mar_KE = Oct2Mar.Q_ne(indxLat_KE,indxLon_KE);
% Q_ne_Apr2Sep_KE = Apr2Sep.Q_ne(indxLat_KE,indxLon_KE);

Q_ce_Oct2Mar_GS = Oct2Mar.Q_ce(indxLat_GS,indxLon_GS);
Q_ce_Apr2Sep_GS = Apr2Sep.Q_ce(indxLat_GS,indxLon_GS);
Q_te_Oct2Mar_GS = Oct2Mar.Q_te(indxLat_GS,indxLon_GS);
Q_te_Apr2Sep_GS = Apr2Sep.Q_te(indxLat_GS,indxLon_GS);
% Q_le_Oct2Mar_GS = Oct2Mar.Q_le(indxLat_GS,indxLon_GS); 
% Q_le_Apr2Sep_GS = Apr2Sep.Q_le(indxLat_GS,indxLon_GS); 
% Q_ne_Oct2Mar_GS = Oct2Mar.Q_ne(indxLat_GS,indxLon_GS);
% Q_ne_Apr2Sep_GS = Apr2Sep.Q_ne(indxLat_GS,indxLon_GS);

Q_ce_Oct2Mar_AR = Oct2Mar.Q_ce(indxLat_AR,indxLon_AR);
Q_ce_Apr2Sep_AR = Apr2Sep.Q_ce(indxLat_AR,indxLon_AR);
Q_te_Oct2Mar_AR = Oct2Mar.Q_te(indxLat_AR,indxLon_AR);
Q_te_Apr2Sep_AR = Apr2Sep.Q_te(indxLat_AR,indxLon_AR);
% Q_le_Oct2Mar_AR = Oct2Mar.Q_le(indxLat_AR,indxLon_AR); 
% Q_le_Apr2Sep_AR = Apr2Sep.Q_le(indxLat_AR,indxLon_AR); 
% Q_ne_Oct2Mar_AR = Oct2Mar.Q_ne(indxLat_AR,indxLon_AR);
% Q_ne_Apr2Sep_AR = Apr2Sep.Q_ne(indxLat_AR,indxLon_AR);

Q_ce_Oct2Mar_BC = Oct2Mar.Q_ce(indxLat_BC,indxLon_BC);
Q_ce_Apr2Sep_BC = Apr2Sep.Q_ce(indxLat_BC,indxLon_BC);
Q_te_Oct2Mar_BC = Oct2Mar.Q_te(indxLat_BC,indxLon_BC);
Q_te_Apr2Sep_BC = Apr2Sep.Q_te(indxLat_BC,indxLon_BC);
% Q_le_Oct2Mar_BC = Oct2Mar.Q_le(indxLat_BC,indxLon_BC); 
% Q_le_Apr2Sep_BC = Apr2Sep.Q_le(indxLat_BC,indxLon_BC); 
% Q_ne_Oct2Mar_BC = Oct2Mar.Q_ne(indxLat_BC,indxLon_BC);
% Q_ne_Apr2Sep_BC = Apr2Sep.Q_ne(indxLat_BC,indxLon_BC);

clear Apr2Sep Oct2Mar

% --- KE ---
[N_Q_ce_Apr2Sep_KE,C_Q_ce_Apr2Sep_KE] = hist(Q_ce_Apr2Sep_KE(:),QbinsCenter); 
P_Q_ce_Apr2Sep_KE = N_Q_ce_Apr2Sep_KE./sum(N_Q_ce_Apr2Sep_KE);

[N_Q_ce_Oct2Mar_KE,C_Q_ce_Oct2Mar_KE] = hist(Q_ce_Oct2Mar_KE(:),QbinsCenter); 
P_Q_ce_Oct2Mar_KE = N_Q_ce_Oct2Mar_KE./sum(N_Q_ce_Oct2Mar_KE);

[N_Q_te_Apr2Sep_KE,C_Q_te_Apr2Sep_KE] = hist(Q_te_Apr2Sep_KE(:),QbinsCenter); 
P_Q_te_Apr2Sep_KE = N_Q_te_Apr2Sep_KE./sum(N_Q_te_Apr2Sep_KE);

[N_Q_te_Oct2Mar_KE,C_Q_te_Oct2Mar_KE] = hist(Q_te_Oct2Mar_KE(:),QbinsCenter); 
P_Q_te_Oct2Mar_KE = N_Q_te_Oct2Mar_KE./sum(N_Q_te_Oct2Mar_KE);

% [N_Q_le_Apr2Sep_KE,C_Q_le_Apr2Sep_KE] = hist(Q_le_Apr2Sep_KE(:),QbinsCenter); 
% P_Q_le_Apr2Sep_KE = N_Q_le_Apr2Sep_KE./sum(N_Q_le_Apr2Sep_KE);
% 
% [N_Q_le_Oct2Mar_KE,C_Q_le_Oct2Mar_KE] = hist(Q_le_Oct2Mar_KE(:),QbinsCenter); 
% P_Q_le_Oct2Mar_KE = N_Q_le_Oct2Mar_KE./sum(N_Q_le_Oct2Mar_KE);
% ----------------------------------

% --- GS ---
[N_Q_ce_Apr2Sep_GS,C_Q_ce_Apr2Sep_GS] = hist(Q_ce_Apr2Sep_GS(:),QbinsCenter); 
P_Q_ce_Apr2Sep_GS = N_Q_ce_Apr2Sep_GS./sum(N_Q_ce_Apr2Sep_GS);

[N_Q_ce_Oct2Mar_GS,C_Q_ce_Oct2Mar_GS] = hist(Q_ce_Oct2Mar_GS(:),QbinsCenter); 
P_Q_ce_Oct2Mar_GS = N_Q_ce_Oct2Mar_GS./sum(N_Q_ce_Oct2Mar_GS);

[N_Q_te_Apr2Sep_GS,C_Q_te_Apr2Sep_GS] = hist(Q_te_Apr2Sep_GS(:),QbinsCenter); 
P_Q_te_Apr2Sep_GS = N_Q_te_Apr2Sep_GS./sum(N_Q_te_Apr2Sep_GS);

[N_Q_te_Oct2Mar_GS,C_Q_te_Oct2Mar_GS] = hist(Q_te_Oct2Mar_GS(:),QbinsCenter); 
P_Q_te_Oct2Mar_GS = N_Q_te_Oct2Mar_GS./sum(N_Q_te_Oct2Mar_GS);
% [N_Q_le_Apr2Sep_GS,C_Q_le_Apr2Sep_GS] = hist(Q_le_Apr2Sep_GS(:),QbinsCenter); 
% P_Q_le_Apr2Sep_GS = N_Q_le_Apr2Sep_GS./sum(N_Q_le_Apr2Sep_GS);
% 
% [N_Q_le_Oct2Mar_GS,C_Q_le_Oct2Mar_GS] = hist(Q_le_Oct2Mar_GS(:),QbinsCenter); 
% P_Q_le_Oct2Mar_GS = N_Q_le_Oct2Mar_GS./sum(N_Q_le_Oct2Mar_GS);
% -------------------------

% --- Agulhus Return currents (AR) ---
[N_Q_ce_Apr2Sep_AR,C_Q_ce_Apr2Sep_AR] = hist(Q_ce_Apr2Sep_AR(:),QbinsCenter); 
P_Q_ce_Apr2Sep_AR = N_Q_ce_Apr2Sep_AR./sum(N_Q_ce_Apr2Sep_AR);

[N_Q_ce_Oct2Mar_AR,C_Q_ce_Oct2Mar_AR] = hist(Q_ce_Oct2Mar_AR(:),QbinsCenter); 
P_Q_ce_Oct2Mar_AR = N_Q_ce_Oct2Mar_AR./sum(N_Q_ce_Oct2Mar_AR);

[N_Q_te_Apr2Sep_AR,C_Q_te_Apr2Sep_AR] = hist(Q_te_Apr2Sep_AR(:),QbinsCenter); 
P_Q_te_Apr2Sep_AR = N_Q_te_Apr2Sep_AR./sum(N_Q_te_Apr2Sep_AR);

[N_Q_te_Oct2Mar_AR,C_Q_te_Oct2Mar_AR] = hist(Q_te_Oct2Mar_AR(:),QbinsCenter); 
P_Q_te_Oct2Mar_AR = N_Q_te_Oct2Mar_AR./sum(N_Q_te_Oct2Mar_AR);

% [N_Q_le_Apr2Sep_AR,C_Q_le_Apr2Sep_AR] = hist(Q_le_Apr2Sep_AR(:),QbinsCenter); 
% P_Q_le_Apr2Sep_AR = N_Q_le_Apr2Sep_AR./sum(N_Q_le_Apr2Sep_AR);
% 
% [N_Q_le_Oct2Mar_AR,C_Q_le_Oct2Mar_AR] = hist(Q_le_Oct2Mar_AR(:),QbinsCenter); 
% P_Q_le_Oct2Mar_AR = N_Q_le_Oct2Mar_AR./sum(N_Q_le_Oct2Mar_AR);
% ----------------------------

% --- Brazil Current (BC) ---
[N_Q_ce_Apr2Sep_BC,C_Q_ce_Apr2Sep_BC] = hist(Q_ce_Apr2Sep_BC(:),QbinsCenter); 
P_Q_ce_Apr2Sep_BC = N_Q_ce_Apr2Sep_BC./sum(N_Q_ce_Apr2Sep_BC);

[N_Q_ce_Oct2Mar_BC,C_Q_ce_Oct2Mar_BC] = hist(Q_ce_Oct2Mar_BC(:),QbinsCenter); 
P_Q_ce_Oct2Mar_BC = N_Q_ce_Oct2Mar_BC./sum(N_Q_ce_Oct2Mar_BC);

[N_Q_te_Apr2Sep_BC,C_Q_te_Apr2Sep_BC] = hist(Q_te_Apr2Sep_BC(:),QbinsCenter); 
P_Q_te_Apr2Sep_BC = N_Q_te_Apr2Sep_BC./sum(N_Q_te_Apr2Sep_BC);

[N_Q_te_Oct2Mar_BC,C_Q_te_Oct2Mar_BC] = hist(Q_te_Oct2Mar_BC(:),QbinsCenter); 
P_Q_te_Oct2Mar_BC = N_Q_te_Oct2Mar_BC./sum(N_Q_te_Oct2Mar_BC);

% [N_Q_le_Apr2Sep_BC,C_Q_le_Apr2Sep_BC] = hist(Q_le_Apr2Sep_BC(:),QbinsCenter); 
% P_Q_le_Apr2Sep_BC = N_Q_le_Apr2Sep_BC./sum(N_Q_le_Apr2Sep_BC);
% 
% [N_Q_le_Oct2Mar_BC,C_Q_le_Oct2Mar_BC] = hist(Q_le_Oct2Mar_BC(:),QbinsCenter); 
% P_Q_le_Oct2Mar_BC = N_Q_le_Oct2Mar_BC./sum(N_Q_le_Oct2Mar_BC);
% ----------------------------
% ======================


%% === make pics ===
f1=figure('Renderer','painters'); % Q

  pic1_size = [8.5 5]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
%   set(f1,'units','normalized','position',[0,0,1,1])
  font_size = 10;
  x_limits = [-17 47];
  y_limits = [0 0.8];
  y_ticks =  [0:0.2:0.8];
  x4text = -12;
  y4text = 0.75;
  ylabel_position=[-0.23, 0.5, 0];
  centerQ = C_Q_ce_Apr2Sep_BC;
  
% ~~~ generate subplot position ~~~  
  row_num=2;        col_num=4;
  margin_left=0.09;  margin_right=0.03;
  margin_top=0.05;  margin_botm=0.1;
  pics_dist_x=0.05; pics_dist_y=0.08;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

   
  subplot('Position',sbplt_posit(1,:))
    hb_ce=bar(centerQ,P_Q_ce_Apr2Sep_KE,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Apr2Sep_KE,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    legend([hb_ce,hb_te],{'Q_{total}^{50m}','Q_{TE}^{50m}'}, ...
        'fontsize',font_size-2);legend boxoff;
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);
    title('KE','fontsize',font_size);
    ylabel(['Apr-Sep',sprintf('\n'),'[%]'],'Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(2,:))
    hb_ce=bar(centerQ,P_Q_ce_Apr2Sep_GS,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Apr2Sep_GS,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);   
    title('GS','fontsize',font_size);
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold'); 
    
  subplot('Position',sbplt_posit(3,:))
    hb_ce=bar(centerQ,P_Q_ce_Apr2Sep_AR,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Apr2Sep_AR,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);   
    title('AR','fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold'); 
    
  subplot('Position',sbplt_posit(4,:))
    hb_ce=bar(centerQ,P_Q_ce_Apr2Sep_BC,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Apr2Sep_BC,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);   
    title('BC','fontsize',font_size);
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold'); 
    
  subplot('Position',sbplt_posit(5,:))
    hb_ce=bar(centerQ,P_Q_ce_Oct2Mar_KE,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Oct2Mar_KE,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);
    ylabel(['Oct-Mar',sprintf('\n'),'[%]'],'Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    xlabel('[W/m^2]','fontsize',font_size);
    text(x4text,y4text,'e','fontsize',font_size,'FontWeight','bold');
    
  subplot('Position',sbplt_posit(6,:))
    hb_ce=bar(centerQ,P_Q_ce_Oct2Mar_GS,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Oct2Mar_GS,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);
    text(x4text,y4text,'f','fontsize',font_size,'FontWeight','bold');
    xlabel('[W/m^2]','fontsize',font_size);
    
  subplot('Position',sbplt_posit(7,:))
    hb_ce=bar(centerQ,P_Q_ce_Oct2Mar_AR,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Oct2Mar_AR,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);
    text(x4text,y4text,'g','fontsize',font_size,'FontWeight','bold'); 
    xlabel('[W/m^2]','fontsize',font_size);    
    
  subplot('Position',sbplt_posit(8,:))
    hb_ce=bar(centerQ,P_Q_ce_Oct2Mar_BC,'r','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);hold on;
    hb_te=bar(centerQ,P_Q_te_Oct2Mar_BC,'b','FaceAlpha',0.2, ...
      'EdgeAlpha',.8,'linewidth',1.2);grid on; 
    set(gca,'xlim',x_limits,'ylim',y_limits,'ytick',y_ticks,'fontsize',font_size);
    text(x4text,y4text,'h','fontsize',font_size,'FontWeight','bold'); 
    xlabel('[W/m^2]','fontsize',font_size);
% ====================

  
%% === output data ===  
print(f1,'-dpng',pic1)
% ====================
