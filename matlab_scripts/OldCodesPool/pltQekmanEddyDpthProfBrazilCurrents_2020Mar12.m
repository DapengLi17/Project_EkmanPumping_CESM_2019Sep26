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
   
% infile1 = '../data_after_manipulation/QekmanEddyGlobalBorealSummerAv_2020Feb11.mat';
% infile2 = '../data_after_manipulation/QekmanEddyGlobalBorealWinterAv_2020Feb11.mat';
infile1 = '../data_after_manipulation/QekmanEddyGlobalDpthProfApr-SepAv_2020Feb11.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalDpthProfOct-MarAv_2020Feb11.mat';

outfile1 = ['../data_after_manipulation/QekmanEddyWinterSummerAvGulfStream_tst_',date_str,'.txt'];
pic1 = ['../pics/QekmanEddyWinterSummerAvGulfStream_',date_str,'.png'];

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

   lat_limits_BC = [-55 -30];  
   lon_limits_BC = [290 320];
%=============================


%% === load data ===
DpthProf_su = load(infile1);
DpthProf_wt = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = DpthProf_su.lat_1d;
lon_1d_raw = DpthProf_su.lon_1d;

indxLat_BC = find(lat_1d_raw >= lat_limits_KE(1) & lat_1d_raw <= lat_limits_KE(2));
lat_1d = lat_1d_raw(indxLat_BC);
indxLon_BC = find(lon_1d_raw >= lon_limits_KE(1) & lon_1d_raw <= lon_limits_KE(2));
lon_1d = lon_1d_raw(indxLon_BC);

z = DpthProf_su.z_wvel_ce; 

Q_li_prof_wt_BC = DpthProf_wt.Q_li(indxLon_BC,indxLat_BC,:); 
Q_li_prof_su_BC = DpthProf_su.Q_li(indxLon_BC,indxLat_BC,:); 
Q_nl_prof_wt_BC = DpthProf_wt.Q_nl(indxLon_BC,indxLat_BC,:);
Q_nl_prof_su_BC = DpthProf_su.Q_nl(indxLon_BC,indxLat_BC,:);
Q_ce_prof_wt_BC = DpthProf_wt.Q_ce(indxLon_BC,indxLat_BC,:);
Q_ce_prof_su_BC = DpthProf_su.Q_ce(indxLon_BC,indxLat_BC,:);

clear Wt Su
  
for iz = 1 : length(z)
   iz

   dummy_Q_ce_prof_wt = Q_ce_prof_wt_BC(:,:,iz);
   [Q_ce_prof_wt_min(iz),Q_ce_prof_wt_av(iz),Q_ce_prof_wt_max(iz)] = ...
      bootstrap5(dummy_Q_ce_prof_wt(:));
   clear dummy_Q_ce_prof_wt 

   dummy_Q_ce_prof_su = Q_ce_prof_su_BC(:,:,iz); 
   [Q_ce_prof_su_min(iz),Q_ce_prof_su_av(iz),Q_ce_prof_su_max(iz)] = ...
       bootstrap5(dummy_Q_ce_prof_su(:));
   clear dummy_Q_ce_prof_su

    dummy_Q_nl_prof_wt = Q_nl_prof_wt_BC(:,:,iz);
    [Q_nl_prof_wt_min(iz),Q_nl_prof_wt_av(iz),Q_nl_prof_wt_max(iz)] = ...
       bootstrap5(dummy_Q_nl_prof_wt(:));
    clear dummy_Q_nl_prof_wt 
    dummy_Q_nl_prof_su = Q_nl_prof_su_BC(:,:,iz);

    [Q_nl_prof_su_min(iz),Q_nl_prof_su_av(iz),Q_nl_prof_su_max(iz)] = ...
       bootstrap5(dummy_Q_nl_prof_su(:));
    clear dummy_Q_nl_prof_su 

end
% ======================


%% === make pics ===
f1=figure('Renderer','painters');  
  % set pic size
  pic1_size=[6.62 6]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
  font_size=10;
  ylabel_position=[-0.10, 0.5, 0]; 
  
% ~~~ generate subplot position ~~~   
  row_num=2; col_num=3; margin_left=0.10;
  margin_right=0.03;margin_top=0.05;margin_botm=0.09;
  pics_dist_x=0.06;pics_dist_y=0.07;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);

% figure format  
%          Global-av  KE  GS BC
% winter
% summer
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

subplot('position',sbplt_posit(1,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_wt_min',Q_ce_prof_wt_av', ...
     Q_ce_prof_wt_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_nl_prof_wt_min',Q_nl_prof_wt_av', ...
     Q_nl_prof_wt_max',-z,'b');grid on;
plot(Q_ce_prof_wt_av',-z,'r*-');hold on;
plot(Q_nl_prof_wt_av',-z,'b*-');grid on;

subplot('position',sbplt_posit(4,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_su_min',Q_ce_prof_su_av', ...
     Q_ce_prof_su_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_nl_prof_su_min',Q_nl_prof_su_av', ...
     Q_nl_prof_su_max',-z,'b');grid on;
plot(Q_ce_prof_su_av',-z,'r*-');hold on;
plot(Q_nl_prof_su_av',-z,'b*-');grid on;
% ==================
