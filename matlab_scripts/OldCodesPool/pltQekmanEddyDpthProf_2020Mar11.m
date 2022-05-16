%% === readme ===

% descrip: matlab scripts plot depth profiles of Qeddy 
% and Qekman for western boudnary currents and global averages. 

% update history:
% v1.0 DL 2020Mar10

% extra notes:
% ===============


%% === set up environments ===
clear all;close all;clc;

  date_str = '2020Mar10';

infile1 = '../data_after_manipulation/QekmanEddyGlobalDpthProfOct-MarAv_2020Feb11.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalDpthProfApr-SepAv_2020Feb11.mat';

  pic1 = ['../pics/VerifyQeddyDpthProf_',date_str,'.png'];
  pic2 = ['../pics/VerifyQekmanDpthProf_',date_str,'.png'];
  
addpath(genpath('Func4EkmanProject/'))

  lon_limits = [0 360];
  lat_limits = [-70 70];
% =============================


%% === load data ===
DpthProf_wt = load(infile1);
DpthProf_su = load(infile2);
% ==================


%% === data analysis ===
lon_1d = DpthProf_wt.lon_1d;
lat_1d = DpthProf_wt.lat_1d;
z = DpthProf_wt.z_wvel_ce; 
% z = = 0 10 20 30 40 50 70 90 110
Q_ce_prof_wt = DpthProf_wt.Q_ce;
Q_ce_prof_su = DpthProf_su.Q_ce;

Q_nl_prof_wt = DpthProf_wt.Q_nl;
Q_nl_prof_su = DpthProf_su.Q_nl;

disp(['working on Dpth prof, large data set, be patient'])
for iz = 1 : length(z)
   iz

   dummy_Q_ce_prof_wt = Q_ce_prof_wt(:,:,iz);
   [Q_ce_prof_wt_min(iz),Q_ce_prof_wt_av(iz),Q_ce_prof_wt_max(iz)] = ...
      bootstrap5(dummy_Q_ce_prof_wt(:));
   clear dummy_Q_ce_prof_wt 

   dummy_Q_ce_prof_su = Q_ce_prof_su(:,:,iz); 
   [Q_ce_prof_su_min(iz),Q_ce_prof_su_av(iz),Q_ce_prof_su_max(iz)] = ...
       bootstrap5(dummy_Q_ce_prof_su(:));
   clear dummy_Q_ce_prof_su

    dummy_Q_nl_prof_wt = Q_nl_prof_wt(:,:,iz);
    [Q_nl_prof_wt_min(iz),Q_nl_prof_wt_av(iz),Q_nl_prof_wt_max(iz)] = ...
       bootstrap5(dummy_Q_nl_prof_wt(:));
    clear dummy_Q_nl_prof_wt 
    dummy_Q_nl_prof_su = Q_nl_prof_su(:,:,iz);

    [Q_nl_prof_su_min(iz),Q_nl_prof_su_av(iz),Q_nl_prof_su_max(iz)] = ...
       bootstrap5(dummy_Q_nl_prof_su(:));
    clear dummy_Q_nl_prof_su 

end
% =========================


%% === make pics ===
f1=figure('Renderer','painters');  
  % set pic size
  pic1_size=[6.62 6]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
  font_size=10;
  ylabel_position=[-0.10, 0.5, 0]; 
  
% ~~~ generate subplot position ~~~   
  row_num=2; col_num=4; margin_left=0.10;
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

subplot('position',sbplt_posit(2,:));
  [h_fil h_plt_Np]=plt_shade_line(Q_ce_prof_su_min',Q_ce_prof_su_av', ...
     Q_ce_prof_su_max',-z,'r');hold on;
  [h_fil h_plt_Np]=plt_shade_line(Q_nl_prof_su_min',Q_nl_prof_su_av', ...
     Q_nl_prof_su_max',-z,'b');grid on;
plot(Q_ce_prof_su_av',-z,'r*-');hold on;
plot(Q_nl_prof_su_av',-z,'b*-');grid on;
% ==================
