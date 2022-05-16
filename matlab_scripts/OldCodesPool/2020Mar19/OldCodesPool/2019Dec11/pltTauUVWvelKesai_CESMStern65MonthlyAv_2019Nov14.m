%% ===readme===

% descrip: matlab scripts plot tau, UVW vel, Kesai from 
% 1) CESM model output
% 2) Stern nonlinear Ekman vertical pumping velocity
% 3) linear Ekman vertical pumping velocity

% update history:
% v1.0 DL 2019Nov14

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Nov14';
  
infile1 = '../data_after_manipulation/CESMtauuvwtempKE_2019Nov06.mat';

addpath(genpath('Func4EkmanProject/'))
rho0 = 1020;
%=============================


%% === load data ===
KE = load(infile1); % KE : Kuroshio Extension
%===================


%% === data analysis ===
% convert lat and lon to x and y
Lat_vec=KE.Lat_vec;Lon_vec=KE.Lon_vec;
Lon=KE.Lon;Lat=KE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);

u = KE.u;
v = KE.v;
w = KE.w;
taux = KE.taux;
tauy = KE.tauy;
tau_A = sqrt(taux.^2+tauy.^2);
tauA_av = squeeze(nanmean(nanmean(tauA,1),2));
jultime_vec = KE.jultime_vec;
jultime= datenum(jultime_vec);
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
for i = 1 : size(u,3) % loop through time
    kesai(:,:,i) = CalcRelVorticityFunc(x,y,u(:,:,i),v(:,:,i));
end

% compute Rossby number
for i = 1 : size(kesai,3)
    Ro(:,:,i) = kesai(:,:,i)./f;  
end

% compute nonlinear Ekman Wvel (Stern 1965)
for i = 1 : size(kesai,3)
 [w_nl(:,:,i)] = CalcEkmanWvelStern65Func(rho0,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,kesai(:,:,i));
end

% compute linear Ekman Wvel
for i = 1 : size(kesai,3)
 [w_li(:,:,i)] = CalcEkmanWvelStern65Func(rho0,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,0);
end

w_ed = w_nl - w_li; % ed: eddy, eddy induced vel
  
% === time average ===
  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% --- monthly av ---
for i = 1 : 12 

%      % taux
%    tauxMonAv(:,:,i) = Calc3dArrayTimeAvFunc(taux,IndxMon{i},'nanmean');     
%      % tauy
%    tauyMonAv(:,:,i) = Calc3dArrayTimeAvFunc(tauy,IndxMon{i},'nanmean'); 

   tau_ampMonAv(:,:,i) = Calc3dArrayTimeAvFunc(tau_A,IndxMon{i},'nanmean'); 
     % u vel
   uMonAv(:,:,i) = Calc3dArrayTimeAvFunc(u,IndxMon{i},'nanmean');    
     % v vel
   vMonAv(:,:,i) = Calc3dArrayTimeAvFunc(v,IndxMon{i},'nanmean');     
     % relative vorticity
   kesaiMonAv(:,:,i) = Calc3dArrayTimeAvFunc(kesai,IndxMon{i},'nanmean');   
     % Rossby num
   RoMonAv(:,:,i) = Calc3dArrayTimeAvFunc(Ro,IndxMon{i},'nanmean');    
     % linear Ekman Wvel
   w_liMonAv(:,:,i) = Calc3dArrayTimeAvFunc(w_li,IndxMon{i},'bootstrap');
     % nonlinear Ekman wvel
   w_nlMonAv(:,:,i) = Calc3dArrayTimeAvFunc(w_nl,IndxMon{i},'bootstrap');
     % eddy W vel (nonlinear-linear Ekman Wvel)
   w_edMonAv(:,:,i) = Calc3dArrayTimeAvFunc(w_ed,IndxMon{i},'bootstrap');
     % CESM W vel
   w_ceMonAv(:,:,i) = Calc3dArrayTimeAvFunc(w,IndxMon{i},'bootstrap');

end

% --- seasonal av ---
% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{12}]; % Winter
  IndxSea{2} = [IndxMon{3};IndxMon{4};IndxMon{5}]; % Spring 
  IndxSea{3} = [IndxMon{6};IndxMon{7};IndxMon{8}]; % Summer
  IndxSea{4} = [IndxMon{9};IndxMon{10};IndxMon{11}]; % Fall
  
for i = 1 : 4 

%      % taux
%    tauxSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(taux,IndxSea{i},'nanmean');     
%      % tauy
%    tauySeaAv(:,:,i) = Calc3dArrayTimeAvFunc(tauy,IndxSea{i},'nanmean');

   tau_ampSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(tau_A,IndxSea{i},'nanmean');
     % u vel
   uSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(u,IndxSea{i},'nanmean');    
     % v vel
   vSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(v,IndxSea{i},'nanmean');     
     % relative vorticity
   kesaiSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(kesai,IndxSea{i},'nanmean');   
     % Rossby num
   RoSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(Ro,IndxSea{i},'nanmean');    
     % linear Ekman Wvel
   w_liSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(w_li,IndxSea{i},'bootstrap');
     % nonlinear Ekman wvel
   w_nlSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(w_nl,IndxSea{i},'bootstrap');
     % eddy W vel (nonlinear-linear Ekman Wvel)
   w_edSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(w_ed,IndxSea{i},'bootstrap');
     % CESM W vel
   w_ceSeaAv(:,:,i) = Calc3dArrayTimeAvFunc(w,IndxSea{i},'bootstrap');

end  
% ------------
% ==================    
  

%% === make pics ===  

  W_lim = [-1e-5 1e-5];
  indxX = 1 : 5 : length(Lon_vec);
  indxY = 1 : 5 : length(Lat_vec);

% monthly av plots  
for i = 1 : size(w_liMonAv,3) 
    
 f1=figure; % winter
  set(f1,'units','normalized','position',[0,0,1,1])
    
 subplot(4,2,2)
   quiver(Lon(indxY,indxX),Lat(indxY,indxX), ...
      uMonAv(indxY,indxX,i),vMonAv(indxY,indxX,i), ...
      'b','AutoScale','off');hold on;
   quiver(132,44,1,0,'b','AutoScale','off');
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
   text(134,44,'1 m/s');
   title(['CESM UV vel [m/s] dpth:',num2str(KE.dpth_u),'m']);
  
 subplot(4,2,3)
  pcolor(Lon,Lat,kesaiMonAv(:,:,i));shading interp;
  caxis([-5 5].*1e-5);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['relative vorticity']);
   
 subplot(4,2,4)
  pcolor(Lon,Lat,RoMonAv(:,:,i));shading interp;
  caxis([-0.5 0.5]);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['Rossby Num: kesai/f']);
  
 subplot(4,2,5)
  pcolor(Lon,Lat,w_liMonAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(linear) [m/s]']);
  
 subplot(4,2,6)
  pcolor(Lon,Lat,w_nlMonAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear) [m/s]']);
  
 subplot(4,2,7)
  pcolor(Lon,Lat,w_edMonAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear)-W(linear) [m/s]']);
  
 subplot(4,2,8)
  pcolor(Lon,Lat,w_ceMonAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);
  
 ax1=subplot(4,2,1);
   pcolor(Lon,Lat,tau_ampMonAv(:,:,i));shading interp;
   caxis([0 0.25]);colorbar;colormap(ax1,parula);
   set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
   title(['CESM wind stress [N/m^2], Mon:',num2str(i,'%02d')])    
  
 pic = ['../pics/TauUVWvelKesai_CESMStern65_Mon',num2str(i,'%02d'),'_',date_str,'.png'];
 print(f1,'-dpng',pic)
 
end

close all

% season av plots  
for i = 1 : size(w_liSeaAv,3) 
    
 f1=figure; % winter
  set(f1,'units','normalized','position',[0,0,1,1])
    
 subplot(4,2,2)
   quiver(Lon(indxY,indxX),Lat(indxY,indxX), ...
      uSeaAv(indxY,indxX,i),vSeaAv(indxY,indxX,i), ...
      'b','AutoScale','off');hold on;
   quiver(132,44,1,0,'b','AutoScale','off');
   set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
   text(134,44,'1 m/s');
   title(['CESM UV vel [m/s] dpth:',num2str(KE.dpth_u),'m']);
  
 subplot(4,2,3)
  pcolor(Lon,Lat,kesaiSeaAv(:,:,i));shading interp;
  caxis([-5 5].*1e-5);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['relative vorticity']);
   
 subplot(4,2,4)
  pcolor(Lon,Lat,RoSeaAv(:,:,i));shading interp;
  caxis([-0.5 0.5]);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['Rossby Num: kesai/f']);
  
 subplot(4,2,5)
  pcolor(Lon,Lat,w_liSeaAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(linear) [m/s]']);
  
 subplot(4,2,6)
  pcolor(Lon,Lat,w_nlSeaAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear) [m/s]']);
  
 subplot(4,2,7)
  pcolor(Lon,Lat,w_edSeaAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar;
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(nonlinear)-W(linear) [m/s]']);
  
 subplot(4,2,8)
  pcolor(Lon,Lat,w_ceSeaAv(:,:,i));shading interp;
  caxis(W_lim);polarmap(64,1);colorbar; 
  set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
  title(['W(cesm) [m/s] dpth:',num2str(KE.dpth_w),'m']);
  
 ax1=subplot(4,2,1);
   pcolor(Lon,Lat,tau_ampSeaAv(:,:,i));shading interp;
   caxis([0 0.25]);colorbar;colormap(ax1,parula);
   set(gca,'xlim',[130.5 169.5],'ylim',[26.5 45.5]);
   title(['CESM wind stress [N/m^2], Sea:',num2str(i,'%02d')])    
  
 pic = ['../pics/TauUVWvelKesai_CESMStern65_Season',num2str(i,'%02d'),'_',date_str,'.png'];
 print(f1,'-dpng',pic)
 
end
% ===================