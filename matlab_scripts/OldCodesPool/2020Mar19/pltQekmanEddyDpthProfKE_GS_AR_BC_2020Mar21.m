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
   
infile1 = '../data_after_manipulation/QekmanEddyDpthProfGlobalAprSepAv_2020Mar19.mat';
infile2 = '../data_after_manipulation/QekmanEddyDpthProfGlobalOctMarAv_2020Mar19.mat';

pic1 = ['../pics/QekmanEddyApr2Sep_Oct2MarAvDpthProfKE_GS_AR_BC_',date_str,'.png'];
pic2 = ['../pics/QekmanApr2Sep_Oct2MarAvDpthProfKE_GS_AR_BC_',date_str,'.png'];

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

% --- Gulf stream region ---
%    lat_limits_GS = [25 50];  
%    lon_limits_GS = [275 325];
% --- good choice ---
%    lat_limits_GS = [32 42];  
%    lon_limits_GS = [285 310];
% -------------------  
   lat_limits_GS = [35 42];  
   lon_limits_GS = [285 305];
% ---------------------------

% --- Kuroshio extension region ---
   lat_limits_KE = [30 45]; 
   lon_limits_KE = [140 170];
% ----------------------------------

% --- Agulhas return currents ---
   lat_limits_AR = [-45 -35];   
   lon_limits_AR = [10 60];
% -----------------------------

% --- Brazil Current region ---
   lat_limits_BC = [-50 -35];  
   lon_limits_BC = [302 315];
% -----------------------------

%=============================


%% === load data ===
DpthProf_Apr2Sep = load(infile1);
DpthProf_Oct2Mar = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = DpthProf_Apr2Sep.lat_1d;
lon_1d_raw = DpthProf_Apr2Sep.lon_1d;

z = DpthProf_Apr2Sep.z_wvel_ce; 

indxLat_GS = find(lat_1d_raw >= lat_limits_GS(1) & lat_1d_raw <= lat_limits_GS(2));
indxLon_GS = find(lon_1d_raw >= lon_limits_GS(1) & lon_1d_raw <= lon_limits_GS(2));

indxLat_KE = find(lat_1d_raw >= lat_limits_KE(1) & lat_1d_raw <= lat_limits_KE(2));
indxLon_KE = find(lon_1d_raw >= lon_limits_KE(1) & lon_1d_raw <= lon_limits_KE(2));

indxLat_BC = find(lat_1d_raw >= lat_limits_BC(1) & lat_1d_raw <= lat_limits_BC(2));
indxLon_BC = find(lon_1d_raw >= lon_limits_BC(1) & lon_1d_raw <= lon_limits_BC(2));

indxLat_AR = find(lat_1d_raw >= lat_limits_AR(1) & lat_1d_raw <= lat_limits_AR(2));
indxLon_AR = find(lon_1d_raw >= lon_limits_AR(1) & lon_1d_raw <= lon_limits_AR(2));

Q_ce_prof_Oct2Mar_GS = DpthProf_Oct2Mar.Q_ce(indxLon_GS,indxLat_GS,:);
Q_ce_prof_Apr2Sep_GS = DpthProf_Apr2Sep.Q_ce(indxLon_GS,indxLat_GS,:);
Q_te_prof_Oct2Mar_GS = DpthProf_Oct2Mar.Q_te(indxLon_GS,indxLat_GS,:);
Q_te_prof_Apr2Sep_GS = DpthProf_Apr2Sep.Q_te(indxLon_GS,indxLat_GS,:);
Q_le_prof_Oct2Mar_GS = DpthProf_Oct2Mar.Q_le(indxLon_GS,indxLat_GS,:); 
Q_le_prof_Apr2Sep_GS = DpthProf_Apr2Sep.Q_le(indxLon_GS,indxLat_GS,:); 
% Q_ne_prof_Oct2Mar_GS = DpthProf_Oct2Mar.Q_ne(indxLon_GS,indxLat_GS,:);
% Q_ne_prof_Apr2Sep_GS = DpthProf_Apr2Sep.Q_ne(indxLon_GS,indxLat_GS,:);
Q_ne_prof_Oct2Mar_GS = Q_te_prof_Oct2Mar_GS - Q_le_prof_Oct2Mar_GS;
Q_ne_prof_Apr2Sep_GS = Q_te_prof_Apr2Sep_GS - Q_le_prof_Apr2Sep_GS;


Q_ce_prof_Oct2Mar_KE = DpthProf_Oct2Mar.Q_ce(indxLon_KE,indxLat_KE,:);
Q_ce_prof_Apr2Sep_KE = DpthProf_Apr2Sep.Q_ce(indxLon_KE,indxLat_KE,:);
Q_te_prof_Oct2Mar_KE = DpthProf_Oct2Mar.Q_te(indxLon_KE,indxLat_KE,:);
Q_te_prof_Apr2Sep_KE = DpthProf_Apr2Sep.Q_te(indxLon_KE,indxLat_KE,:);
Q_le_prof_Oct2Mar_KE = DpthProf_Oct2Mar.Q_le(indxLon_KE,indxLat_KE,:); 
Q_le_prof_Apr2Sep_KE = DpthProf_Apr2Sep.Q_le(indxLon_KE,indxLat_KE,:); 
% Q_ne_prof_Oct2Mar_KE = DpthProf_Oct2Mar.Q_ne(indxLon_KE,indxLat_KE,:);
% Q_ne_prof_Apr2Sep_KE = DpthProf_Apr2Sep.Q_ne(indxLon_KE,indxLat_KE,:);
Q_ne_prof_Oct2Mar_KE = Q_te_prof_Oct2Mar_KE - Q_le_prof_Oct2Mar_KE;
Q_ne_prof_Apr2Sep_KE = Q_te_prof_Apr2Sep_KE - Q_le_prof_Apr2Sep_KE;

  
Q_ce_prof_Oct2Mar_BC = DpthProf_Oct2Mar.Q_ce(indxLon_BC,indxLat_BC,:);
Q_ce_prof_Apr2Sep_BC = DpthProf_Apr2Sep.Q_ce(indxLon_BC,indxLat_BC,:);
Q_te_prof_Oct2Mar_BC = DpthProf_Oct2Mar.Q_te(indxLon_BC,indxLat_BC,:);
Q_te_prof_Apr2Sep_BC = DpthProf_Apr2Sep.Q_te(indxLon_BC,indxLat_BC,:);
Q_le_prof_Oct2Mar_BC = DpthProf_Oct2Mar.Q_le(indxLon_BC,indxLat_BC,:); 
Q_le_prof_Apr2Sep_BC = DpthProf_Apr2Sep.Q_le(indxLon_BC,indxLat_BC,:); 
% Q_ne_prof_Oct2Mar_BC = DpthProf_Oct2Mar.Q_ne(indxLon_BC,indxLat_BC,:);
% Q_ne_prof_Apr2Sep_BC = DpthProf_Apr2Sep.Q_ne(indxLon_BC,indxLat_BC,:);
Q_ne_prof_Oct2Mar_BC = Q_te_prof_Oct2Mar_BC - Q_le_prof_Oct2Mar_BC;
Q_ne_prof_Apr2Sep_BC = Q_te_prof_Apr2Sep_BC - Q_le_prof_Apr2Sep_BC;


Q_ce_prof_Oct2Mar_AR = DpthProf_Oct2Mar.Q_ce(indxLon_AR,indxLat_AR,:);
Q_ce_prof_Apr2Sep_AR = DpthProf_Apr2Sep.Q_ce(indxLon_AR,indxLat_AR,:);
Q_te_prof_Oct2Mar_AR = DpthProf_Oct2Mar.Q_te(indxLon_AR,indxLat_AR,:);
Q_te_prof_Apr2Sep_AR = DpthProf_Apr2Sep.Q_te(indxLon_AR,indxLat_AR,:);
Q_le_prof_Oct2Mar_AR = DpthProf_Oct2Mar.Q_le(indxLon_AR,indxLat_AR,:); 
Q_le_prof_Apr2Sep_AR = DpthProf_Apr2Sep.Q_le(indxLon_AR,indxLat_AR,:); 
% Q_ne_prof_Oct2Mar_AR = DpthProf_Oct2Mar.Q_ne(indxLon_AR,indxLat_AR,:);
% Q_ne_prof_Apr2Sep_AR = DpthProf_Apr2Sep.Q_ne(indxLon_AR,indxLat_AR,:);
Q_ne_prof_Oct2Mar_AR = Q_te_prof_Oct2Mar_AR - Q_le_prof_Oct2Mar_AR;
Q_ne_prof_Apr2Sep_AR = Q_te_prof_Apr2Sep_AR - Q_le_prof_Apr2Sep_AR;

clear DpthProf_su DpthProf_wt


for iz = 1 : length(z)
   iz

% --- Kuroshio extension ---
   dummy_Q_ce_prof_Oct2Mar_KE = Q_ce_prof_Oct2Mar_KE(:,:,iz);
   [Q_ce_prof_Oct2Mar_KE_min(iz),Q_ce_prof_Oct2Mar_KE_av(iz), ...
      Q_ce_prof_Oct2Mar_KE_max(iz)] = bootstrap5(dummy_Q_ce_prof_Oct2Mar_KE(:));

   dummy_Q_ce_prof_Apr2Sep_KE = Q_ce_prof_Apr2Sep_KE(:,:,iz); 
   [Q_ce_prof_Apr2Sep_KE_min(iz),Q_ce_prof_Apr2Sep_KE_av(iz), ...
      Q_ce_prof_Apr2Sep_KE_max(iz)] = bootstrap5(dummy_Q_ce_prof_Apr2Sep_KE(:));
   clear dummy_Q_ce_prof_Oct2Mar_KE dummy_Q_ce_prof_Apr2Sep_KE

   dummy_Q_te_prof_Oct2Mar_KE = Q_te_prof_Oct2Mar_KE(:,:,iz);
    [Q_te_prof_Oct2Mar_KE_min(iz),Q_te_prof_Oct2Mar_KE_av(iz), ...
      Q_te_prof_Oct2Mar_KE_max(iz)] = bootstrap5(dummy_Q_te_prof_Oct2Mar_KE(:));
    
    dummy_Q_te_prof_Apr2Sep_KE = Q_te_prof_Apr2Sep_KE(:,:,iz);
    [Q_te_prof_Apr2Sep_KE_min(iz),Q_te_prof_Apr2Sep_KE_av(iz), ...
      Q_te_prof_Apr2Sep_KE_max(iz)] = bootstrap5(dummy_Q_te_prof_Apr2Sep_KE(:));   
   clear dummy_Q_te_prof_Oct2Mar_KE dummy_Q_te_prof_Apr2Sep_KE
   
    dummy_Q_ne_prof_Oct2Mar_KE = Q_ne_prof_Oct2Mar_KE(:,:,iz);
    [Q_ne_prof_Oct2Mar_KE_min(iz),Q_ne_prof_Oct2Mar_KE_av(iz), ...
      Q_ne_prof_Oct2Mar_KE_max(iz)] = bootstrap5(dummy_Q_ne_prof_Oct2Mar_KE(:));
    
    dummy_Q_ne_prof_Apr2Sep_KE = Q_ne_prof_Apr2Sep_KE(:,:,iz);
    [Q_ne_prof_Apr2Sep_KE_min(iz),Q_ne_prof_Apr2Sep_KE_av(iz), ...
      Q_ne_prof_Apr2Sep_KE_max(iz)] = bootstrap5(dummy_Q_ne_prof_Apr2Sep_KE(:));
   clear dummy_Q_ne_prof_Oct2Mar_KE dummy_Q_ne_prof_Apr2Sep_KE
   
   dummy_Q_le_prof_Oct2Mar_KE = Q_le_prof_Oct2Mar_KE(:,:,iz);
    [Q_le_prof_Oct2Mar_KE_min(iz),Q_le_prof_Oct2Mar_KE_av(iz), ...
      Q_le_prof_Oct2Mar_KE_max(iz)] = bootstrap5(dummy_Q_le_prof_Oct2Mar_KE(:));
    
    dummy_Q_le_prof_Apr2Sep_KE = Q_le_prof_Apr2Sep_KE(:,:,iz);
    [Q_le_prof_Apr2Sep_KE_min(iz),Q_le_prof_Apr2Sep_KE_av(iz), ...
      Q_le_prof_Apr2Sep_KE_max(iz)] = bootstrap5(dummy_Q_le_prof_Apr2Sep_KE(:));
   clear dummy_Q_le_prof_Oct2Mar_KE dummy_Q_le_prof_Apr2Sep_KE
% ---------------------

% --- Gulf Stream ---
   dummy_Q_ce_prof_Oct2Mar_GS = Q_ce_prof_Oct2Mar_GS(:,:,iz);
   [Q_ce_prof_Oct2Mar_GS_min(iz),Q_ce_prof_Oct2Mar_GS_av(iz), ...
      Q_ce_prof_Oct2Mar_GS_max(iz)] = bootstrap5(dummy_Q_ce_prof_Oct2Mar_GS(:));

   dummy_Q_ce_prof_Apr2Sep_GS = Q_ce_prof_Apr2Sep_GS(:,:,iz); 
   [Q_ce_prof_Apr2Sep_GS_min(iz),Q_ce_prof_Apr2Sep_GS_av(iz), ...
      Q_ce_prof_Apr2Sep_GS_max(iz)] = bootstrap5(dummy_Q_ce_prof_Apr2Sep_GS(:));
   clear dummy_Q_ce_prof_Oct2Mar_GS dummy_Q_ce_prof_Apr2Sep_GS

   dummy_Q_te_prof_Oct2Mar_GS = Q_te_prof_Oct2Mar_GS(:,:,iz);
    [Q_te_prof_Oct2Mar_GS_min(iz),Q_te_prof_Oct2Mar_GS_av(iz), ...
      Q_te_prof_Oct2Mar_GS_max(iz)] = bootstrap5(dummy_Q_te_prof_Oct2Mar_GS(:));
    
    dummy_Q_te_prof_Apr2Sep_GS = Q_te_prof_Apr2Sep_GS(:,:,iz);
    [Q_te_prof_Apr2Sep_GS_min(iz),Q_te_prof_Apr2Sep_GS_av(iz), ...
      Q_te_prof_Apr2Sep_GS_max(iz)] = bootstrap5(dummy_Q_te_prof_Apr2Sep_GS(:));   
   clear dummy_Q_te_prof_Oct2Mar_GS dummy_Q_te_prof_Apr2Sep_GS
   
    dummy_Q_ne_prof_Oct2Mar_GS = Q_ne_prof_Oct2Mar_GS(:,:,iz);
    [Q_ne_prof_Oct2Mar_GS_min(iz),Q_ne_prof_Oct2Mar_GS_av(iz), ...
      Q_ne_prof_Oct2Mar_GS_max(iz)] = bootstrap5(dummy_Q_ne_prof_Oct2Mar_GS(:));
    
    dummy_Q_ne_prof_Apr2Sep_GS = Q_ne_prof_Apr2Sep_GS(:,:,iz);
    [Q_ne_prof_Apr2Sep_GS_min(iz),Q_ne_prof_Apr2Sep_GS_av(iz), ...
      Q_ne_prof_Apr2Sep_GS_max(iz)] = bootstrap5(dummy_Q_ne_prof_Apr2Sep_GS(:));
   clear dummy_Q_ne_prof_Oct2Mar_GS dummy_Q_ne_prof_Apr2Sep_GS
   
   dummy_Q_le_prof_Oct2Mar_GS = Q_le_prof_Oct2Mar_GS(:,:,iz);
    [Q_le_prof_Oct2Mar_GS_min(iz),Q_le_prof_Oct2Mar_GS_av(iz), ...
      Q_le_prof_Oct2Mar_GS_max(iz)] = bootstrap5(dummy_Q_le_prof_Oct2Mar_GS(:));
    
    dummy_Q_le_prof_Apr2Sep_GS = Q_le_prof_Apr2Sep_GS(:,:,iz);
    [Q_le_prof_Apr2Sep_GS_min(iz),Q_le_prof_Apr2Sep_GS_av(iz), ...
      Q_le_prof_Apr2Sep_GS_max(iz)] = bootstrap5(dummy_Q_le_prof_Apr2Sep_GS(:));
   clear dummy_Q_le_prof_Oct2Mar_GS dummy_Q_le_prof_Apr2Sep_GS
% ---------------------

% --- Agulhas Return Current region ---
   dummy_Q_ce_prof_Oct2Mar_AR = Q_ce_prof_Oct2Mar_AR(:,:,iz);
   [Q_ce_prof_Oct2Mar_AR_min(iz),Q_ce_prof_Oct2Mar_AR_av(iz), ...
      Q_ce_prof_Oct2Mar_AR_max(iz)] = bootstrap5(dummy_Q_ce_prof_Oct2Mar_AR(:));

   dummy_Q_ce_prof_Apr2Sep_AR = Q_ce_prof_Apr2Sep_AR(:,:,iz); 
   [Q_ce_prof_Apr2Sep_AR_min(iz),Q_ce_prof_Apr2Sep_AR_av(iz), ...
      Q_ce_prof_Apr2Sep_AR_max(iz)] = bootstrap5(dummy_Q_ce_prof_Apr2Sep_AR(:));
   clear dummy_Q_ce_prof_Oct2Mar_AR dummy_Q_ce_prof_Apr2Sep_AR

   dummy_Q_te_prof_Oct2Mar_AR = Q_te_prof_Oct2Mar_AR(:,:,iz);
    [Q_te_prof_Oct2Mar_AR_min(iz),Q_te_prof_Oct2Mar_AR_av(iz), ...
      Q_te_prof_Oct2Mar_AR_max(iz)] = bootstrap5(dummy_Q_te_prof_Oct2Mar_AR(:));
    
    dummy_Q_te_prof_Apr2Sep_AR = Q_te_prof_Apr2Sep_AR(:,:,iz);
    [Q_te_prof_Apr2Sep_AR_min(iz),Q_te_prof_Apr2Sep_AR_av(iz), ...
      Q_te_prof_Apr2Sep_AR_max(iz)] = bootstrap5(dummy_Q_te_prof_Apr2Sep_AR(:));   
   clear dummy_Q_te_prof_Oct2Mar_AR dummy_Q_te_prof_Apr2Sep_AR
   
    dummy_Q_ne_prof_Oct2Mar_AR = Q_ne_prof_Oct2Mar_AR(:,:,iz);
    [Q_ne_prof_Oct2Mar_AR_min(iz),Q_ne_prof_Oct2Mar_AR_av(iz), ...
      Q_ne_prof_Oct2Mar_AR_max(iz)] = bootstrap5(dummy_Q_ne_prof_Oct2Mar_AR(:));
    
    dummy_Q_ne_prof_Apr2Sep_AR = Q_ne_prof_Apr2Sep_AR(:,:,iz);
    [Q_ne_prof_Apr2Sep_AR_min(iz),Q_ne_prof_Apr2Sep_AR_av(iz), ...
      Q_ne_prof_Apr2Sep_AR_max(iz)] = bootstrap5(dummy_Q_ne_prof_Apr2Sep_AR(:));
   clear dummy_Q_ne_prof_Oct2Mar_AR dummy_Q_ne_prof_Apr2Sep_AR
   
   dummy_Q_le_prof_Oct2Mar_AR = Q_le_prof_Oct2Mar_AR(:,:,iz);
    [Q_le_prof_Oct2Mar_AR_min(iz),Q_le_prof_Oct2Mar_AR_av(iz), ...
      Q_le_prof_Oct2Mar_AR_max(iz)] = bootstrap5(dummy_Q_le_prof_Oct2Mar_AR(:));
    
    dummy_Q_le_prof_Apr2Sep_AR = Q_le_prof_Apr2Sep_AR(:,:,iz);
    [Q_le_prof_Apr2Sep_AR_min(iz),Q_le_prof_Apr2Sep_AR_av(iz), ...
      Q_le_prof_Apr2Sep_AR_max(iz)] = bootstrap5(dummy_Q_le_prof_Apr2Sep_AR(:));
   clear dummy_Q_le_prof_Oct2Mar_AR dummy_Q_le_prof_Apr2Sep_AR
% -----------------

% --- Brazil currents ---
   dummy_Q_ce_prof_Oct2Mar_BC = Q_ce_prof_Oct2Mar_BC(:,:,iz);
   [Q_ce_prof_Oct2Mar_BC_min(iz),Q_ce_prof_Oct2Mar_BC_av(iz), ...
      Q_ce_prof_Oct2Mar_BC_max(iz)] = bootstrap5(dummy_Q_ce_prof_Oct2Mar_BC(:));

   dummy_Q_ce_prof_Apr2Sep_BC = Q_ce_prof_Apr2Sep_BC(:,:,iz); 
   [Q_ce_prof_Apr2Sep_BC_min(iz),Q_ce_prof_Apr2Sep_BC_av(iz), ...
      Q_ce_prof_Apr2Sep_BC_max(iz)] = bootstrap5(dummy_Q_ce_prof_Apr2Sep_BC(:));
   clear dummy_Q_ce_prof_Oct2Mar_BC dummy_Q_ce_prof_Apr2Sep_BC

   dummy_Q_te_prof_Oct2Mar_BC = Q_te_prof_Oct2Mar_BC(:,:,iz);
    [Q_te_prof_Oct2Mar_BC_min(iz),Q_te_prof_Oct2Mar_BC_av(iz), ...
      Q_te_prof_Oct2Mar_BC_max(iz)] = bootstrap5(dummy_Q_te_prof_Oct2Mar_BC(:));
    
    dummy_Q_te_prof_Apr2Sep_BC = Q_te_prof_Apr2Sep_BC(:,:,iz);
    [Q_te_prof_Apr2Sep_BC_min(iz),Q_te_prof_Apr2Sep_BC_av(iz), ...
      Q_te_prof_Apr2Sep_BC_max(iz)] = bootstrap5(dummy_Q_te_prof_Apr2Sep_BC(:));   
   clear dummy_Q_te_prof_Oct2Mar_BC dummy_Q_te_prof_Apr2Sep_BC
   
    dummy_Q_ne_prof_Oct2Mar_BC = Q_ne_prof_Oct2Mar_BC(:,:,iz);
    [Q_ne_prof_Oct2Mar_BC_min(iz),Q_ne_prof_Oct2Mar_BC_av(iz), ...
      Q_ne_prof_Oct2Mar_BC_max(iz)] = bootstrap5(dummy_Q_ne_prof_Oct2Mar_BC(:));
    
    dummy_Q_ne_prof_Apr2Sep_BC = Q_ne_prof_Apr2Sep_BC(:,:,iz);
    [Q_ne_prof_Apr2Sep_BC_min(iz),Q_ne_prof_Apr2Sep_BC_av(iz), ...
      Q_ne_prof_Apr2Sep_BC_max(iz)] = bootstrap5(dummy_Q_ne_prof_Apr2Sep_BC(:));
   clear dummy_Q_ne_prof_Oct2Mar_BC dummy_Q_ne_prof_Apr2Sep_BC
   
   dummy_Q_le_prof_Oct2Mar_BC = Q_le_prof_Oct2Mar_BC(:,:,iz);
    [Q_le_prof_Oct2Mar_BC_min(iz),Q_le_prof_Oct2Mar_BC_av(iz), ...
      Q_le_prof_Oct2Mar_BC_max(iz)] = bootstrap5(dummy_Q_le_prof_Oct2Mar_BC(:));
    
    dummy_Q_le_prof_Apr2Sep_BC = Q_le_prof_Apr2Sep_BC(:,:,iz);
    [Q_le_prof_Apr2Sep_BC_min(iz),Q_le_prof_Apr2Sep_BC_av(iz), ...
      Q_le_prof_Apr2Sep_BC_max(iz)] = bootstrap5(dummy_Q_le_prof_Apr2Sep_BC(:));
   clear dummy_Q_le_prof_Oct2Mar_BC dummy_Q_le_prof_Apr2Sep_BC
% ---------------------


end
% ======================


%% === make pics ===
f1=figure('Renderer','painters');  
  % set pic size
  pic1_size=[6.62 6]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
  font_size=10;
  ylabel_position=[-0.10, 0.5, 0]; 
  x_lim = [0 100];
  x_ticks = [0:20:60];
% ~~~ generate subplot position ~~~   
  row_num=2; col_num=4; margin_left=0.10;
  margin_right=0.03;margin_top=0.05;margin_botm=0.09;
  pics_dist_x=0.06;pics_dist_y=0.07;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);

% figure format  
%           KE  GS  AR  BC
% Apr-Sep
% Oct-Mar
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subplot('position',sbplt_posit(1,:));
  [h_fil h_plt_ce]=plt_shade_line(Q_ce_prof_Apr2Sep_KE_min',Q_ce_prof_Apr2Sep_KE_av', ...
     Q_ce_prof_Apr2Sep_KE_max',-z,'r');hold on;
  [h_fil h_plt_te]=plt_shade_line(Q_te_prof_Apr2Sep_KE_min',Q_te_prof_Apr2Sep_KE_av', ...
     Q_te_prof_Apr2Sep_KE_max',-z,'b');grid on;
  legend([h_plt_te,h_plt_ce],{'Q_{TE}','Q_{total}'});legend boxoff;
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('KE');
ylabel('Apr-Sep');

subplot('position',sbplt_posit(2,:));
  [h_fil h_plt_ce]=plt_shade_line(Q_ce_prof_Apr2Sep_GS_min',Q_ce_prof_Apr2Sep_GS_av', ...
     Q_ce_prof_Apr2Sep_GS_max',-z,'r');hold on;
  [h_fil h_plt_te]=plt_shade_line(Q_te_prof_Apr2Sep_GS_min',Q_te_prof_Apr2Sep_GS_av', ...
     Q_te_prof_Apr2Sep_GS_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('GS');

subplot('position',sbplt_posit(3,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_Apr2Sep_AR_min',Q_ce_prof_Apr2Sep_AR_av', ...
     Q_ce_prof_Apr2Sep_AR_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Apr2Sep_AR_min',Q_te_prof_Apr2Sep_AR_av', ...
     Q_te_prof_Apr2Sep_AR_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('AR');

subplot('position',sbplt_posit(4,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_Apr2Sep_BC_min',Q_ce_prof_Apr2Sep_BC_av', ...
     Q_ce_prof_Apr2Sep_BC_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Apr2Sep_BC_min',Q_te_prof_Apr2Sep_BC_av', ...
     Q_te_prof_Apr2Sep_BC_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('BC');

subplot('position',sbplt_posit(5,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_Oct2Mar_KE_min',Q_ce_prof_Oct2Mar_KE_av', ...
     Q_ce_prof_Oct2Mar_KE_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_KE_min',Q_te_prof_Oct2Mar_KE_av', ...
     Q_te_prof_Oct2Mar_KE_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
ylabel('Oct-Mar');xlabel('[W/m^2]');

subplot('position',sbplt_posit(6,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_Oct2Mar_GS_min',Q_ce_prof_Oct2Mar_GS_av', ...
     Q_ce_prof_Oct2Mar_GS_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_GS_min',Q_te_prof_Oct2Mar_GS_av', ...
     Q_te_prof_Oct2Mar_GS_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
xlabel('[W/m^2]');

subplot('position',sbplt_posit(7,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_Oct2Mar_AR_min',Q_ce_prof_Oct2Mar_AR_av', ...
     Q_ce_prof_Oct2Mar_AR_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_AR_min',Q_te_prof_Oct2Mar_AR_av', ...
     Q_te_prof_Oct2Mar_AR_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
xlabel('[W/m^2]');

subplot('position',sbplt_posit(8,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_Oct2Mar_BC_min',Q_ce_prof_Oct2Mar_BC_av', ...
     Q_ce_prof_Oct2Mar_BC_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_BC_min',Q_te_prof_Oct2Mar_BC_av', ...
     Q_te_prof_Oct2Mar_BC_max',-z,'b');grid on;
set(gca,'xlim',x_lim,'XTick',x_ticks)
xlabel('[W/m^2]');
% ==================


%% === make pics ===
f2=figure('Renderer','painters');  
  % set pic size
  pic1_size=[6.62 6]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f2,'Units','inches','Position',[5,5,pic1_size]);
  font_size=10;
  ylabel_position=[-0.10, 0.5, 0]; 
  x_lim = [-5 10];
  x_ticks = [-5:5:10];
% ~~~ generate subplot position ~~~   
  row_num=2; col_num=4; margin_left=0.10;
  margin_right=0.03;margin_top=0.05;margin_botm=0.09;
  pics_dist_x=0.06;pics_dist_y=0.07;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);

% figure format  
%           KE  GS  AR  BC
% Apr-Sep
% Oct-Mar
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
subplot('position',sbplt_posit(1,:));
  [h_fil h_plt_ce]=plt_shade_line(Q_te_prof_Apr2Sep_KE_min',Q_te_prof_Apr2Sep_KE_av', ...
     Q_te_prof_Apr2Sep_KE_max',-z,'r');hold on;
  [h_fil h_plt_te]=plt_shade_line(Q_le_prof_Apr2Sep_KE_min',Q_le_prof_Apr2Sep_KE_av', ...
     Q_le_prof_Apr2Sep_KE_max',-z,'b');grid on;
  [h_fil h_plt_te]=plt_shade_line(Q_ne_prof_Apr2Sep_KE_min',Q_ne_prof_Apr2Sep_KE_av', ...
     Q_ne_prof_Apr2Sep_KE_max',-z,'k');
  legend([h_plt_te,h_plt_ce],{'Q_{TE}','Q_{total}'});legend boxoff;
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('KE');
ylabel('Apr-Sep');

subplot('position',sbplt_posit(2,:));
  [h_fil h_plt_ce]=plt_shade_line(Q_te_prof_Apr2Sep_GS_min',Q_te_prof_Apr2Sep_GS_av', ...
     Q_te_prof_Apr2Sep_GS_max',-z,'r');hold on;
  [h_fil h_plt_te]=plt_shade_line(Q_le_prof_Apr2Sep_GS_min',Q_le_prof_Apr2Sep_GS_av', ...
     Q_le_prof_Apr2Sep_GS_max',-z,'b');grid on;
  [h_fil h_plt_te]=plt_shade_line(Q_ne_prof_Apr2Sep_GS_min',Q_ne_prof_Apr2Sep_GS_av', ...
     Q_ne_prof_Apr2Sep_GS_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('GS');

subplot('position',sbplt_posit(3,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Apr2Sep_AR_min',Q_te_prof_Apr2Sep_AR_av', ...
     Q_te_prof_Apr2Sep_AR_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_le_prof_Apr2Sep_AR_min',Q_le_prof_Apr2Sep_AR_av', ...
     Q_le_prof_Apr2Sep_AR_max',-z,'b');grid on;
  [h_fil h_plt_Np]=plt_shade_line(Q_ne_prof_Apr2Sep_AR_min',Q_ne_prof_Apr2Sep_AR_av', ...
     Q_ne_prof_Apr2Sep_AR_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('AR');

subplot('position',sbplt_posit(4,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Apr2Sep_BC_min',Q_te_prof_Apr2Sep_BC_av', ...
     Q_te_prof_Apr2Sep_BC_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_le_prof_Apr2Sep_BC_min',Q_le_prof_Apr2Sep_BC_av', ...
     Q_le_prof_Apr2Sep_BC_max',-z,'b');grid on;
  [h_fil h_plt_Np]=plt_shade_line(Q_ne_prof_Apr2Sep_BC_min',Q_ne_prof_Apr2Sep_BC_av', ...
     Q_ne_prof_Apr2Sep_BC_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
title('BC');

subplot('position',sbplt_posit(5,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_KE_min',Q_te_prof_Oct2Mar_KE_av', ...
     Q_te_prof_Oct2Mar_KE_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_le_prof_Oct2Mar_KE_min',Q_le_prof_Oct2Mar_KE_av', ...
     Q_le_prof_Oct2Mar_KE_max',-z,'b');grid on;
  [h_fil h_plt_Np]=plt_shade_line(Q_ne_prof_Oct2Mar_KE_min',Q_ne_prof_Oct2Mar_KE_av', ...
     Q_ne_prof_Oct2Mar_KE_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
ylabel('Oct-Mar');xlabel('[W/m^2]');

subplot('position',sbplt_posit(6,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_GS_min',Q_te_prof_Oct2Mar_GS_av', ...
     Q_te_prof_Oct2Mar_GS_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_le_prof_Oct2Mar_GS_min',Q_le_prof_Oct2Mar_GS_av', ...
     Q_le_prof_Oct2Mar_GS_max',-z,'b');grid on;
  [h_fil h_plt_Np]=plt_shade_line(Q_ne_prof_Oct2Mar_GS_min',Q_ne_prof_Oct2Mar_GS_av', ...
     Q_ne_prof_Oct2Mar_GS_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
xlabel('[W/m^2]');

subplot('position',sbplt_posit(7,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_AR_min',Q_te_prof_Oct2Mar_AR_av', ...
     Q_te_prof_Oct2Mar_AR_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_le_prof_Oct2Mar_AR_min',Q_le_prof_Oct2Mar_AR_av', ...
     Q_le_prof_Oct2Mar_AR_max',-z,'b');grid on;
  [h_fil h_plt_Np]=plt_shade_line(Q_ne_prof_Oct2Mar_AR_min',Q_ne_prof_Oct2Mar_AR_av', ...
     Q_ne_prof_Oct2Mar_AR_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
xlabel('[W/m^2]');

subplot('position',sbplt_posit(8,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_te_prof_Oct2Mar_BC_min',Q_te_prof_Oct2Mar_BC_av', ...
     Q_te_prof_Oct2Mar_BC_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_le_prof_Oct2Mar_BC_min',Q_le_prof_Oct2Mar_BC_av', ...
     Q_le_prof_Oct2Mar_BC_max',-z,'b');grid on;
  [h_fil h_plt_Np]=plt_shade_line(Q_ne_prof_Oct2Mar_BC_min',Q_ne_prof_Oct2Mar_BC_av', ...
     Q_ne_prof_Oct2Mar_BC_max',-z,'k');
set(gca,'xlim',x_lim,'XTick',x_ticks)
xlabel('[W/m^2]');
% ==================

  
%% === output data ===  
print(f1,'-dpng',pic1)
print(f2,'-dpng',pic2)
% ====================
