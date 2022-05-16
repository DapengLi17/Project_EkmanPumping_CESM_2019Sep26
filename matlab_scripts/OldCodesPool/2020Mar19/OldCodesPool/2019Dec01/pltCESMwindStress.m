%% ===readme===

% descrip: matlab scripts extract parameters from CESM nc files, 
% regrid the data and save them to .mat file   

% update history:
% v1.0 DL 2019Oct05

% extra notes:
% parameters needed: u,v,w velocity (UVEL,VVEL,WVEL), 
% wind stress (TAUX, TAUY) and temperature (TEMP). 
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Oct15';
  
%InDir = '/atlantis3/zhangqiuying/regrid_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02/ocn_kuroshio/';
InDir1 = '../raw_data/';
  TAUX_Str = {'TAUX','time','ULONG','ULAT'};% TAUX files
  TAUY_Str = {'TAUY','time','ULONG','ULAT'}; % TAUY files

% Infile = [InDir1 'kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_WVEL_00870101-00871231.nc'];
% ncdisp(Infile)
%============================


%% === load data ===
[jultime,ulat,ulon,taux_raw] = LoadCESMKEncFiles3DVarFunc(InDir1,TAUX_Str);
[~,~,~,tauy_raw] = LoadCESMKEncFiles3DVarFunc(InDir1,TAUY_Str);
%===================

%% ==
% unit of wind stress taux and tauy are dyne/centimeter²
% 1 dyne/centimeter² [dyn/cm²] = 0.1 newton/meter² [N/m²]
% see https://www.translatorscafe.com/unit-converter/en/pressure/25-18/dyne%2Fcentimeter%C2%B2-newton%2Fmeter%C2%B2/
taux_unitconvt = taux_raw.*0.1;
tauy_unitconvt = tauy_raw.*0.1;

for i = 1 : size(taux_unitconvt,3);
    
    clear dummy_taux dummy_tauy dummy_tau
    dummy_taux = taux_unitconvt(:,:,i);
    dummy_tauy = tauy_unitconvt(:,:,i);
    dummy_tau = sqrt(dummy_taux.^2+dummy_tauy.^2);
    tauav(i)=nanmean(dummy_tau(:));
    taumax(i)=nanmax(dummy_tau(:));
    
end

figure;
 plot([1:365],tauav,'b*-');grid on;
 
 indx= find(tauav>0.3)
  tauav(indx)
  tauav(51) 