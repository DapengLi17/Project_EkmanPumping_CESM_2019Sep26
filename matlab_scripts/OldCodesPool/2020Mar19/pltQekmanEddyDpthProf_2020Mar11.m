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

for iz = 1 : length(z)
    iz
    dummy_Q_ce_prof_wt = Q_ce_prof_wt(:,:,iz);
    dummy_Q_ce_prof_su = Q_ce_prof_su(:,:,iz); 
    dummy_Q_nl_prof_wt = Q_nl_prof_wt(:,:,iz);
    dummy_Q_nl_prof_su = Q_nl_prof_su(:,:,iz);
    
  [Q_ce_prof_wt_min(iz),Q_ce_prof_wt_av(iz),Q_ce_prof_wt_max(iz)] = ...
      bootstrap5(dummy_Q_ce_prof_wt(:));
  [Q_ce_prof_su_min(iz),Q_ce_prof_su_av(iz),Q_ce_prof_su_max(iz)] = ...
      bootstrap5(dummy_Q_ce_prof_su(:));
  [Q_nl_prof_wt_min(iz),Q_nl_prof_wt_av(iz),Q_nl_prof_wt_max(iz)] = ...
      bootstrap5(dummy_Q_nl_prof_wt(:));
  [Q_nl_prof_su_min(iz),Q_nl_prof_su_av(iz),Q_nl_prof_su_max(iz)] = ...
      bootstrap5(dummy_Q_nl_prof_su(:));
  
end
