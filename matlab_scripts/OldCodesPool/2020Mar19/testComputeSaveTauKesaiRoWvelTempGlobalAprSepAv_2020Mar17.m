clear all;close all;clc;

infile1 = '../data_after_manipulation/WvelTempGlobal_CESM_87-04_old.mat';
infile2 = '../data_after_manipulation/WvelTempGlobal_CESM_87-04.mat';
pic1 = '../pics/testComputeSaveTauKesaiRoWvelTempGlobalAprSepAv_2020Mar17.png';

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

%%
dummy1 = load(infile1);
dummy2 = load(infile2);

%%
f1=figure;

 set(f1,'units','normalized','position',[0,0,1,1]);
 iday = 1;

m_proj('miller','lon',[0 360],'lat',[-70 70]); 
 
 subplot(2,1,1)
  m_pcolor(dummy1.lon_1d,dummy1.lat_1d,dummy1.w_nl_3d(:,:,iday).*86400);
   shading interp;caxis([-2 2]);colorbar;polarmap;
  m_coast('patch',[.7 .7 .7],'edgecolor','none');
  m_grid('linestyle','none');axis normal;
  title('W_Ekman from ComputeSaveTauKesaiRoWvelTempGlobalBorealSummer_2020Feb08.m [m/day]');
  
 subplot(2,1,2)
  m_pcolor(dummy2.lon_1d,dummy2.lat_1d,dummy2.w_te_3d(:,:,iday).*86400);
   shading interp;caxis([-2 2]);colorbar;polarmap;   
  m_coast('patch',[.7 .7 .7],'edgecolor','none');
  m_grid('linestyle','none');axis normal; 
  title('W_Ekman from ComputeSaveTauKesaiRoWvelTempGlobalAprSepAv_2020Mar17.m [m/day]')
   
print(f1,'-dpng',pic1)

