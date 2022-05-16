function  [pcx, pcy, pr, parea, aerr] = func_fit_circle_area ( xvec, yvec )
%
%    Description : Fit circle to the given closed contour/polygon data and 
%                     do a shape test. Circle fitting based on closed contour
%                     area and centroid, as adopted from Wikipedia.
%
%    Input : 1.  xvec - (vector, float) : x-vertices of polygon/closed contour
%            2.  yvec - (vector, float) : y-vertices of polygon/closed contour
%
%    NOTES:  1. xvec and yvec should same units and same size or length, and
%                  should be for a closed contour/polygon.
%            2. If only one Nx2 variable is passed for arg1 and arg2, 1st 
%                  column is taken as xvec and second column as yvec.  
%            4. the "shape" test will work well only with km-based axes, &
%                  not with lat-lon based axes.
%
%
%    Output : 1. pcx - (scalar, float) : center X of fitted circle
%             2. pcy - (scalar, float) : center Y of fitted circle
%             3. pr  - (scalar, float) : radius of fitted circle
%             4. parea- (scalar, float) : area of the closed contour/circle
%             4. aerr- (scalar, float) : area error for fitted circle
%                                (percentage error in circle's area with respect
%                                  to the "deformed" area of polygon (best if 
%                                  0-40%). This value may depend on value of Q 
%                                  or W parameter, and level of smoothing applied.
%                               Please note that the value for aerr can be higher
%                               100, though it is expressed as %.
%
%            Units same as that of input variables.
%
%            All output variables are set to NaN, if any one of the following
%            is true :
%                 1. number of contour data points < 3
%                 2. bad values in contour data
%                 3. failure in the area (or shape) test
%   
%
%   Jaison Kurian (jaison@atmos.ucla.edu)
%   Aug/26/2009
%

%-----------
% Initialize
%-----------

    ndims = sum(size(xvec) > 1) ;      % for compatibility with 
    if ( nargin == 1 && ndims == 2 )   %   original CIRCLEFIT style of
       tmp  = xvec ;                   %   input  
       xvec = tmp(:,1);                %   
       yvec = tmp(:,2);
    elseif ( nargin < 2 )
       error ([mfilename ':argchk'], '\n\t FATAL (%s) : Need 2 args input, if input vars are vectors.\n',mfilename)
    end
    npts = numel(xvec) ;
    ccx = NaN ; ccy = NaN ; r = NaN ; aerr = NaN ; % initialize output variables to NaN
    
    if ( npts < 3 )  % need at least 3-points to fit a circle    
       return        % output variables = NaN
    end 
    
    % use dbmod to activate/deactivate this part
    
    if ( numel(yvec) ~= npts )
       error ([mfilename ':argchk'], '\n\t FATAL (%s) : Size/Length mismatch between X-points (%03.0d) \& Y-points (%03.0d).\n',mfilename, npts, numel(yvec))
    end
    
    % dbmod region over
    
    xvec = xvec(:) ; % force column vector
    yvec = yvec(:) ; %         "

%-----------
% Shape Test   : Read all 
%-----------
%   Method : 
%
%            Area = 1/2 * sum_i=0_to_n-1_of ( (x_i * y_i+1) - (x_i+1 * y_i) )
%
%            Centroid_x = 1/(6*area) * sum_i=0_to_n-1_of ( (x_i + x_i+1) * 
%                                 ( (x_i * y_i+1) - (x_i+1 * y_i) ) )
%            Centroid_y = 1/(6*area) * sum_i=0_to_n-1_of ( (y_i + y_i+1) * 
%                                 ( (x_i * y_i+1) - (x_i+1 * y_i) ) )
%
%   NOTE : For practical purpose, "polygon" implies closed contour.
%
%   Ref : http://en.wikipedia.org/wiki/Polygon
%

    % area and centroid of closed contour/polygon

    pii   = 4 * atan (1) ;   

    tmp   =  ( xvec(1:npts-1) .* yvec(2:npts) ) - ( xvec(2:npts) .* yvec(1:npts-1) ) ;

    polAr = 1/2 * sum ( tmp );
    parea = abs ( polAr )    ;  % force +ve values
    pr    = sqrt (parea/pii) ; % radius from area

    pcx = 1/(6*polAr) * sum ( ( xvec(1:npts-1) + xvec(2:npts) ) .* tmp ) ;
    pcy = 1/(6*polAr) * sum ( ( yvec(1:npts-1) + yvec(2:npts) ) .* tmp ) ;

    % Centroid test here is skipped, since the "area" test is robust to
    %    locate filament structures

    %dist_cent = ( ( pcx - ccx )^2 + ( pcy - ccy)^2 )^0.5 ;  
    %cent_err  = 100.0 * dist_cent/r ; 

    %if ( cent_err > mperr_cent )          % if centroid test fails
     %  ccx = NaN ; ccy = NaN ; r = NaN ; aerr = NaN ; % return with NaN values
     %  return
    %end      

    %------------------------------------------------- 
    % area test (to skip non-circular shaped contours)
    %------------------------------------------------- 

          %  Find distance between circle center and contour points
    dist_poly = ( ( xvec - pcx ).^2 + ( yvec - pcy ).^2 ).^0.5 ;

          %  for polygon points outside the circle, take that on the circle
          %
          %              C 
          %             /|
          %            / |       A --> (ccx, ccy)    circle center
          %           /  |       C --> (xvec, yvec)  polygon point
          %         D+   |       D --> Point on the circle, corresponding to C 
          %         /|   |
          %        / |   |       
          %       /  |   |       sin_theta (A) = BC/AC = ED/AD = ED/r
          %      /   |   |       --> ED = r * sin_theta
          %     /    |   |       --> circle_y cy = ccy + ED
          %    +---------+
          %   A      E    B
          %
          %  Note that this procudre will not create "self intersecting" polygons.
          %     if the contour's top and bottom lies outside the circle, then for 
          %     thos points there will be a "line" instead of "negative areas" and
          %     in this situation, the simple poly-area method will work just fine. 
          %     (test with x = [-3 3 2  1 -1 -2 -3]; y =[ 0 0 3  0  0  3  0] ) 
     
    sintheta = (yvec - pcy)./dist_poly ;
    ptmp_y   = pcy + ( pr .* sintheta ) ;
    ptmp_x   = pcx - (pcx - xvec) .* (  (pcy - ptmp_y)./(pcy - yvec) )  ; 

    pout_id  = find (dist_poly > pr) ; % indices of polygon points outside circle
                       % p_inon_? : polygon x or y points inside & on the circle
    p_inon_x = xvec ;                 % init 
    p_inon_y = yvec ;                 % init 
    p_inon_x(pout_id) = ptmp_x(pout_id) ;
    p_inon_y(pout_id) = ptmp_y(pout_id) ;

         % area of closed contour/polygon enclosed by the circle
    tmp   =  ( p_inon_x(1:npts-1) .* p_inon_y(2:npts) ) ...
                   - ( p_inon_x(2:npts) .* p_inon_y(1:npts-1) ) ;
    parea_incirc = 1/2 * abs ( sum (tmp) ) ; 

         % area error 
         %
         %    parea_incirc  - area of polygon enclosed in fitted circle
         %    parea_outcirc - area of polygon outside fitted circle
         %    
         %    if the circle fit is 100% accurate, then the following 
         %    two relations will hold : 
         %
         %            parea_incirc       = carea
         %       -->  parea_incirc/carea = 1
         %       -->  1 - parea_incirc/carea = 0 -----(1)
         %
         %            parea_outcirc       = 0         
         %       -->  parea_outcirc/carea = 0    -----(2) 
         %    
         %    (1)+(2) --> ( 1 - parea_incirc/carea ) + parea_outcirc/carea = 0 --(3)
         %  
         %    If the fitting is not accurate, then (3) will not be zero, and will
         %    give a realistic estimate of the fitting error, normalized by the 
         %    area of fitted circle. It also avoids the need for centroid test.
         % 
         %    Please note that error estimate using (3) can be higher than 100. 
         % 
    aerr = 100.0 * ( ( 1 - parea_incirc/parea ) + ( parea - parea_incirc )/parea ) ; 

    return
