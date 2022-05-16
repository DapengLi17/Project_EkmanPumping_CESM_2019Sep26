Q_nl_df(abs(Q_nl_df)<5)=nan;
Q_ce_df(abs(Q_ce_df)<5)=nan;

ratioQnl2ce=Q_nl_df./Q_ce_df;
ratioQnl2ce(ratioQnl2ce>2) = NaN;
ratioQnl2ce(ratioQnl2ce<-2) = NaN;

figure
 pcolor(Lon,Lat,ratioQnl2ce);shading interp;caxis([-1 1]);colorbar;polarmap
figure
 pcolor(Lon,Lat,Ro_wt_std);shading interp;caxis([-1 1]); 
 
 
 Qlimt=20;
 k1 = find(Q_nl_df > Qlimt);
 k2 = find(Q_ce_df > Qlimt);
 
 k = intersect(k1,k2);length(k)
  
 figure;
  plot(Q_nl_df(k));hold on;plot(Q_ce_df(k))

  [rQ_nl2ce,P] = corrcoef(Q_nl_df(k),Q_ce_df(k))
  


figure
 subplot(4,1,1);
  plot(Q_nl_df_Ro,'b');hold on;plot(Q_nl_df_Ro_Dspk,'r');grid on;
 subplot(4,1,2)
  plot(Q_ce_df_Ro,'b');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;
 subplot(4,1,3)
  plot(Q_nl_df_Ro_Dspk,'b');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;
 subplot(4,1,4)
  plot(Q_nl_df_Ro_Dspk,'b*');hold on;plot(Q_ce_df_Ro_Dspk,'r*');grid on;
  
Q_nl_df_Ro_Dspk(abs(Q_nl_df_Ro_Dspk)<10)=nan;
Q_ce_df_Ro_Dspk(abs(Q_ce_df_Ro_Dspk)<10)=nan;

indxNoNaN_Qnl = find(~isnan(Q_nl_df_Ro_Dspk));
indxNoNaN_Qce = find(~isnan(Q_ce_df_Ro_Dspk));
indxNoNaN_Q = intersect(indxNoNaN_Qnl,indxNoNaN_Qce); length(indxNoNaN_Q)

[rQ_nl2ce,P] = corrcoef(Q_nl_df_Ro_Dspk(indxNoNaN_Q),Q_ce_df_Ro_Dspk(indxNoNaN_Q))


