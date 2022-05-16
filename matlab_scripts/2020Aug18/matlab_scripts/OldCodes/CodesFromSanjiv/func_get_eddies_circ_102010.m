function [eddy emask] = func_get_eddies_circ(etype, grid, diag, cont_vals, min_radius, max_radius, mperr_area, min_points,min_slaamp, emask)
%
% Description : Get the following details of closed contours of specified value :
%                   - center position (x0,y0)
%                   - radius
%                   - vorticity sign
%                   - mean vorticity, temp, salt, spiciness, sigma_theta
%                   - mean 
%
% Input   : etype   - string, scalar, eddy type ('cyc' or 'anti')
%           grid    - struct, created by calling func_init_grid
%           diag    - struct, created by calling func_get_diag
%           cont_vals- vector, values of contour to search for
%           min_radius - minimum radius required to define an eddy
%           max_radius - maximum radius required to define an eddy
%           mperr_area - maximum percentage arror in area (circle fit)
%           min_points - required number of minimum grid points inside a 
%                             closed contour
%           min_slaamp  - float, scalar, minimum SLA (cm) amplitude
%           emask      - mask with NaN values over those eddis which are
%                             identified
%            
%
% Output  : cyc.x  - cyclone, center x position
%           cyc.y  - cyclone, center y position
%           cyc.r  - cyclone, center y position
%           ncycl  - number of cyclones 
%              -----> all above fields for anticyclone
%
% Written By : Jaison Kurian (jaison@atmos.ucla.edu)
%              Xavier Capet  (xavier@atmos.ucla.edy)
%              Inputs from Francois Colas (francois@atmos.ucla.edu)
%
% Written On : May/24/2008
%
%  Modifications : July/30/2008 : Basic changes to suit new func_fit_circle and
%                                   its "shape" test.
%                                 Irrespective of the area test or minimum radius,
%                                   the circle with respect to vorticity should 
%                                   have a minimum of 5 points inside it.
%                                 Capability of handling snapshots without eddies.
%                                 Minimum radius is applied to radius based on vort.
%                  Sep/18/2008 : Work based on diag.sla now. The eddy "identification"
%                                 from SLA maps with ellipse fitting have significant
%                                 input from Xavier Capet. The idea of normalizing SLA
%                                 has been suggested by Francois Colas.
%                                 No serious stuff like radius estimation based on vorticity or 
%                                eddy interior averaging.
%                  Sep/18/2008 : Ellipse fitting instead of circle.
%
%
%                  June/01/2009 : - need a minimum number of 8 grid points inside the closed
%                          contour (in addition to the minimum radius of 50 km)
%                    - Absolute value of SLA within the closed contour (even on
%                         single grid point) shouldn't be lower than that at the 
%                         boundary in any case. That is, SLA values increase 
%                         monotonically towards the eddy  center. (sla_test)
%
%                    - All "in on" test is done for the closed contour instead of 
%                          the fitted circle.
%
%         June/22/2009 - New eddy identification procedure :
%                           - closed contours can have both +ve and -ve values 
%                               inside it, now the sign test is using gradient
%                           - one type of eddy at a time, with contour values from
%                                -ve extreme to +ve extreme for anticyclones and
%                                exactly the opposite way for cyclones. They need
%                                to be supplied in that shape.
%                           - got rid of sla_test and introduced eddy amplitude
%                                 test (minimum 1 cm)
%         Aug/19/2009 - added a foolproof eddy sign test
%                        removed global variable vnames, and introduced
%                          vdef (variables_default).
%
%        Aug/19/2009 : radius is from area enclosed by the contour, center is 
%                           the center of a circle fitted to the closed contour
%
%        Oct/04/2010 : func_fit_circle is replaced by func_fit_circle_area. The
%                       only difference between these functions is "trimming"
%                       unwanted lines and there is abolutely no difference in
%                       terms of results. See modification above Aug/19/2009.
%                       - func_poly_area is not used anymore
%                       - input argument for min_slaamp
%------------------------------------------------------------------------

% initialize

     if ( nargin < 9 )
          error([mfilename ':argchk'], '\n   FATAL (%s) : Too few input arguments. \n', mfilename)
     end

     etype = lower(etype) ;
     switch etype
        case 'cyc'
           cont_vals = sort([cont_vals],2,'descend') ;
        case 'anti'
           cont_vals = sort([cont_vals],2,'ascend') ;
        otherwise
           error('\n\t FATAL (%s) : Eddy type (arg 1) should be either cyc or anti.\n', mfilename)
     end

     ncvals   = length (cont_vals) ;
     neddy    = 0 ;
     epsilon  = 1E-24  ;               % check to avoid divisions by zero
     eddy     = [] ; 
     pii     = 4*atan(1) ;

     % NOTE order of contours : anticyclones, -ve extreme to +ve extreme
     %                          cyclones      +ve extreme to -ve extreme   

      
     % NOTE : successive contours : For the very first contour value, eddy 
     %           identification depends only on closed contours and eddy tests.
     %           From the second contour onwards, additional test is made to 
     %           see whether this eddy is already identified or not.

     if ( nargin == 9 ) ; % emask is not supplied
         emask = diag.sla * 0 + 1;
     end


     for icv = 1:ncvals ;

    
         diag.sla = diag.sla .* emask ;

         cont_val = [cont_vals(icv) cont_vals(icv)];
         cont_c_km   = contour(grid.xkm, grid.ykm, diag.sla, cont_val) ;
         cont_c_ll   = contour(grid.lon, grid.lat, diag.sla, cont_val) ;

%         cont_val
%         'size cont_c_km',size(cont_c_km)
%         'size cont_c_ll',size(cont_c_ll) 
%         cont_c_km(:,1:end)
%         cont_c_ll(:,1:end)

%         pause 

         if ( size (cont_c_km,2) > 3 ) 


            cont_d_km   = func_get_contour_data(cont_c_km);
            cont_d_ll   = func_get_contour_data(cont_c_ll);
            ncont       = size(cont_d_km);          % ncont - number of contours 
            ncont       = ncont(2) ; 

%            cont_d_km
%            cont_d_ll 
%            ncont
%            pause 

             
            for ico  = 1:ncont ;              % ico - iterate through contours       

              diag.sla = diag.sla .* emask ; % need to update for possible inner contours 
              if ( cont_d_km(ico).isopen == 0 && numel(cont_d_km(ico).xdata) > 3 ) 

                  cont_x_km  = cont_d_km(ico).xdata;  % contour xkm points
                  cont_y_km  = cont_d_km(ico).ydata;  % contour ykm points
                  cont_x_ll  = cont_d_ll(ico).xdata;  % contour lon points
                  cont_y_ll  = cont_d_ll(ico).ydata;  % contour lat points
                 
                  [cent_x_km, cent_y_km, radius_km, Carea, errkm] = func_fit_circle_area(cont_x_km,cont_y_km) ; % center & radius of contour, in km
                  [cent_x_ll, cent_y_ll, radius_ll, Carea_ll, errll] = func_fit_circle_area(cont_x_ll,cont_y_ll) ; % center & radius of contour, in lat-lon
                 
                  idkm = [] ; coid = [] ;
                  idkm        = func_gridid (grid.xkm, grid.ykm, cent_x_km, cent_y_km) ;
                  yx0km       = [cent_y_km, cent_x_km] ;

                  sla_c         =  diag.sla(idkm(1),idkm(2)) ;
                            
                  % tests : 1-pass and 0-fail
           
                  shape_test  = 0 ;  minrad_test = 0  ; sign_test  = 0 ; ampl_test = 0 ;
                  maxrad_test = 0 ;  ocean_test  = 0  ; minpts_test = 0 ; nan_test = 0 ;
                  
                  if ( errkm <= mperr_area )     ; shape_test = 1 ; end
                  if ( radius_km >= min_radius ) ; minrad_test = 1 ; end
                  if ( radius_km <= max_radius ) ; maxrad_test = 1 ; end
                  if ( ~isnan(sla_c) &&  ~isempty(sla_c)) ; ocean_test = 1; end
                  
                  clear in on inon sla_inon ;
                  [in on]      =  func_inpoly (grid.xkm, grid.ykm, cont_x_km, cont_y_km) ;
                  
                  inon         =  max(in,on) ;
                  sla_inon     = diag.sla(inon) ;
                  sla_nan      = numel(sla_inon(isnan(sla_inon))) ;


                  if ( sla_nan == 0 )
                      nan_test = 1 ;
                  end

                  if ( numel(sla_inon) >= min_points )
                      minpts_test = 1 ; 
                  end 

                  sla_max = NaN ; % this is going to be a tricky definition, watch out...

                  if ( strcmp(etype,'cyc') == 1 )
                      sla_max      =  min(diag.sla(in))      ; % Max is the min for CYC
                      ampl         =  cont_vals(icv)-sla_max ; % amplitude is always +ve 
                      if ( max(diag.sla(in)) > cont_vals(icv) ) ; ampl = -10000 ; end      ;
                  else                                         % for anticyclones
                      sla_max      =  max(diag.sla(in))      ; % Max is the max for ANTI
                      ampl         = sla_max-cont_vals(icv)  ; % amplitude is always +ve
                      if ( min(diag.sla(in)) < cont_vals(icv) ) ; ampl = -10000 ; end      ;
                  end 

                  if ( ampl >= min_slaamp )
                       sign_test = 1 ;
                  end 

%                  [shape_test, minrad_test, maxrad_test, ocean_test, minpts_test, sign_test, nan_test] 
         
                  eddy_test = shape_test * minrad_test * maxrad_test * ocean_test *  minpts_test  *  sign_test  * nan_test;

           
                  if ( eddy_test == 1 ) 

                      clear tmp1 tmp2
          
                      tmp1        =  squeeze ( grid.xkm(idkm(1),:) ) ;
                      tmp2        =  squeeze ( grid.mask(idkm(1),:) ) ;
                      dco_km      =  func_dfcoast ( tmp1, tmp2, idkm(2) ) ;  % distance from coast                      
                      sla_ave     =  tsnanmean(sla_inon) ;
                      spd_max     =  max(diag.spd(inon));
% Added by Sanjiv (6/15/2020)
                      sst_inon    =  diag.sst(inon);
                      sst_ave     =  tsnanmean(sst_inon);
                      sss_inon    =  diag.sss(inon);
                      sss_ave     =  tsnanmean(sss_inon);
                      utau_inon   =  diag.utau(inon);
                      utau_ave    =  tsnanmean(utau_inon);
                      vtau_inon   =  diag.vtau(inon);
                      vtau_ave    =  tsnanmean(vtau_inon);
                      qnet_inon   =  diag.qnet(inon);
                      qnet_ave    =  tsnanmean(qnet_inon);
                      sbl_inon    =  diag.sbl(inon);
                      sbl_ave     =  tsnanmean(sbl_inon);
%%
%%
                      vort_inon   =  diag.vort(inon) ;
                      vort_ave    =  tsnanmean(vort_inon) ;
                      if ( strcmp(etype,'cyc') == 1 )
                          vort_max = max(diag.vort(in)) ; 
                      else
                          vort_max = min(diag.vort(in)) ; 
                      end
                      neddy = neddy + 1 ;
    
                      eddy.r(neddy)     = radius_ll ;
                      eddy.x(neddy)     = cent_x_ll ;
                      eddy.y(neddy)     = cent_y_ll ;
                      eddy.i(neddy)     = idkm(2)   ;
                      eddy.j(neddy)     = idkm(1)   ;
                      eddy.rkm(neddy)   = radius_km ;
                      eddy.xkm(neddy)   = cent_x_km ;
                      eddy.ykm(neddy)   = cent_y_km ;
                      eddy.dcokm(neddy) = dco_km    ;
                      eddy.slamax(neddy)= sla_max   ;
                      eddy.slaave(neddy)= sla_ave   ;
% Added by Sanjiv
                      eddy.sstave(neddy) = sst_ave   ;
                      eddy.sssave(neddy) = sss_ave   ;
                      eddy.utau(  neddy)= utau_ave   ;
                      eddy.vtau(  neddy)= vtau_ave   ;
                      eddy.qnet(  neddy)= qnet_ave   ;
                      eddy.sbl(   neddy)=  sbl_ave   ;
%%
                      eddy.cont(neddy)  = cont_vals(icv) ;
                      eddy.perr(neddy)  = errkm     ;
                      eddy.area(neddy)  = Carea     ;          % 14
                      eddy.vort(neddy)  = vort_ave  ;
                      eddy.vortmax(neddy)= vort_max  ;
                      eddy.swirlmax(neddy)= spd_max  ;
                      eddy.dkm(neddy)   = NaN       ; % supply with nans now. 
                      eddy.d(neddy)     = NaN       ; % supply with nans now. 

                      emask(inon)       = NaN ; % mask the grid points corresponding to identified eddy

                   %   fprintf ('\t  %s  icv = %4i cont = %7.2f  ico = %4i  Npt = %6i \n', etype, icv, cont_vals(icv), ico, numel(cont_x_ll) )
                  end  % if (eddy_test)
              end  % if (cont_d_km(ico).isopen ==0) 
          end  % for ico=1:ncont 
       end  % if (size(cont_c_km,2) > 3) 
     end  % for icv=1:ncvals 

     vdef     = {'r' 'x' 'y' 'i' 'j' 'rkm' 'xkm' 'ykm' 'dcokm' 'area' 'slaave' 'slamax' 'cont' 'perr' 'd' 'dkm'} ;
     % if no eddies are found, then supply NaN values

     nvar   = length(vdef) ;

%     if ( neddy == 0 )
%        neddy  = 1 ;
%        for iv = 1:nvar ;
%           e_mat   = ['eddy.', vdef{iv}] ;
%           eval( [e_mat, ' = NaN',';' ]) ;
%        end
%     end

     eddy.n  = neddy; 

     return

