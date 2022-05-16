pic2 = ['../pics/QContourWinterStorm_KE_',date_str,'.png'];

%% --- Q winter storms ---
f2=figure('Renderer','painters');
 
  % set pic size
  pic2_size=[12 10]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f2,'Units','inches','Position',[5,5,pic2_size]);
  font_size=10;
  ylabel_position=[-0.08, 0.5, 0];
  Xlim     = [130 170];     Xticks=[130:10:170];
  Ylim     = [26.5 45.5];   Yticks=[30:5:45];
  Q_lim_ce = [-400 400];    Q_ce_Ticks = [-400:200:400];
  Q_lim_nl = [-200 200];    Q_nl_Ticks = [-200:100:200];
  
% ~~~ generate subplot position ~~~  
  row_num=3;col_num=3;margin_left=0.06;
  margin_right=0.07;margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.03; pics_dist_y=0.05;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

 subplot('Position',sbplt_posit(1,:))
  pcolor(Lon,Lat,Q_li_ns);shading interp;caxis(Q_lim_nl);
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size);     
    ylabel(['Q_f',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    title(['storm-free days av (N: ',num2str(length(indx_ns) ),' days)']);
    text(132,28,'a','fontsize',font_size,'FontWeight','bold');
        
 subplot('Position',sbplt_posit(2,:))
  pcolor(Lon,Lat,Q_li_st);shading interp;caxis(Q_lim_nl);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'yticklabel',[],'fontsize',font_size);    
    title(['storm days av (N: ',num2str(length(indx_st) ),' days)']);
    text(132,28,'b','fontsize',font_size,'FontWeight','bold');  
    
 subplot('Position',sbplt_posit(3,:))
  pcolor(Lon,Lat,Q_li_st-Q_li_ns);shading interp;caxis(Q_lim_nl);
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);    
    title(['storm days av - storm-free days av']);
    text(132,28,'c','fontsize',font_size,'FontWeight','bold');  
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(3,2) .01 sbplt_posit(3,4)-0.02], ...
       'ylim',Q_lim_nl,'YTick',Q_nl_Ticks);
    title(hc,'[W/m^2]','FontSize',font_size-2);
  
 subplot('Position',sbplt_posit(4,:))
  pcolor(Lon,Lat,Q_nl_ns);shading interp;caxis(Q_lim_nl);
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size);    
    ylabel(['Q_{\xi}',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    text(132,28,'d','fontsize',font_size,'FontWeight','bold'); 
    
 subplot('Position',sbplt_posit(5,:))
  pcolor(Lon,Lat,Q_nl_st);shading interp;caxis(Q_lim_nl);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'yticklabel',[],'fontsize',font_size);    
    text(132,28,'e','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(6,:))
  pcolor(Lon,Lat,Q_nl_st-Q_nl_ns);shading interp;caxis(Q_lim_nl);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'yticklabel',[],'fontsize',font_size);    
    text(132,28,'f','fontsize',font_size,'FontWeight','bold'); 
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',Q_lim_nl,'YTick',Q_nl_Ticks);
    title(hc,'[W/m^2]','FontSize',font_size-2);
    
 subplot('Position',sbplt_posit(7,:))
  pcolor(Lon,Lat,Q_ce_ns);shading interp;caxis(Q_lim_ce);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size); 
    xlabel('Lon','fontsize',font_size);
    ylabel(['Q_{t}',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    text(132,28,'g','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(8,:))
  pcolor(Lon,Lat,Q_ce_st);shading interp;caxis(Q_lim_ce);
     hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);
    xlabel('Lon','fontsize',font_size);
    text(132,28,'h','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(9,:))
  pcolor(Lon,Lat,Q_ce_st-Q_ce_ns);shading interp;caxis(Q_lim_ce);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);   
    set(gca,'xlim',Xlim,'XTick',Xticks,'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size); 
    xlabel('Lon','fontsize',font_size);
    text(132,28,'i','fontsize',font_size,'FontWeight','bold');
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(9,2) .01 sbplt_posit(9,4)-0.02], ...
       'ylim',Q_lim_ce,'YTick',Q_ce_Ticks);
    title(hc,'[W/m^2]','FontSize',font_size-2);
    
  polarmap(64,1);
% ====================


%% === output data ===
print(f2,'-dpng','-r500',pic2)
% ====================
