   figure(fno) ; colormap(blue_red) ;
      set (gcf,'Position', [20 470 620 650]) ;
      [hC hC] = contourf(grid.lon,grid.lat,diag.sla,[-50:5:200])  ; hold on ;
      set (hC, 'LineStyle', 'none') ; axis equal ; axis ([120 200 23 50]) ;
      caxis([-50 200]) ; % very important
      colorbar ;
      
      for ied = 1:antipr.n
        [anti_cx anti_cy]  = func_get_circle (antipr.x(ied), antipr.y(ied), antipr.r(ied)) ;
        plot (anti_cx, anti_cy, '-k','LineWidth', 1.75)
      end
      
      for ied = 1:cycpr.n
        [cyc_cx cyc_cy]  = func_get_circle (cycpr.x(ied), cycpr.y(ied), cycpr.r(ied)) ;
        plot (cyc_cx, cyc_cy, '-r','LineWidth',1.75)
      end
      
      hold off;


