function vonz = func_zsurf(var, zfield, zval, dbm)
%
% Description : Extract variable at specified z-surface (depth or isopycnal). 
%                  This is a slightly modified version of sigma2z.m taken from
%                  Francois Colas (belongs to ROMS TOOLS or not??).
%
% Input  : 1. var    (3D ZYX Matrix) - variable to be extracted 
%          2. zfield (3D ZYX Matrix) - extract based on this field  (can be +ve/-ve)
%          3. zval   (scalar)        - extract var on zfield=zval surface  (    "  )
%          4. dbm    (scalar,optional) - debug mode, 1 for debugging and 0 for not.
%                                       (overrides the global variable "dbmod").
%
% Global  : dbmod  (scalar) - debug mode, 1 for debugging and 0 for not
%                                                                            
% Output : vonz (2D YX Matrix): variable on specified Z-surface
%
%
% Based on   : sigma2z.m ( c/o Francois Colas )
% Written By : UCLA ROMS Team (jaison@atmos.ucla.edu)
% Written On : June/11/2008
% Copyright @ UCLA ROMS Team
% Tool       : Eddy Tracker
% Version    : 1.0
%
%----------------------------------------------------------------
  global dbmod
%-----------
% Initialize
%-----------

  if ( ~exist('dbm', 'var') ) ; dbm = dbmod ; end  % choose debugging mode
 
  if ( dbm == 1 ) % if debugging mode on
      if ( nargin < 3 ) 
          error ([mfilename ':argchk'], '\n\t FATAL (%s) : Need 3 input arguments (var, zfield, zval).\n', mfilename) 
      else   
         %if ( ~func_shape(var,3) || ~func_shape(zfield,3) ||  numel(var) ~= numel(zfield) )
         if ( numel(var) ~= numel(zfield) )
           error ([mfilename ':varchk'], '\n\t FATAL (%s) : Arguments 1 and 2 (var and zfield) should be of same size/shape.\n', mfilename) ;
         end
         if ( ~isscalar(zval) )
           error ([mfilename ':varchk'], '\n\t FATAL (%s) : Argument 3 (zval) should be scalar.\n', mfilename) ;
         end
      end
  end

%--------------
% Preprocessing
%--------------

  zfield      =  squeeze(abs(zfield))     ;
  zval        =  squeeze(abs(zval))       ;
  nz          =  size(zfield,1)           ;
 
  vex         =  zfield > zval   ;     % 1  where z > zval
  levs        =  squeeze(sum(vex,1));  % number of sigma levels where z > depth
  levs(levs==nz) = nz - 1;             % fix for surface values
  warning off ; mask = levs ./ levs ; warning on ; % NaN if levs = 0

  if ( sum( size(var)> 1) == 3 )  % ZYX-data
     [nz,ny,nx]  =  size(zfield)    ;
     [imat,jmat] = meshgrid((1:nx),(1:ny)); % 2D matrix with i and j indices if var
     pos         = nz*ny*(imat-1)+nz*(jmat-1)+levs; % linear position
  elseif ( sum( size(var)> 1) == 2 )
     [nz,nx]     =  size(zfield)    ;
     imat        = [1:1:nx] ; % 1D matrix with i indices of var
     pos         = nz*(imat-1)+levs; % linear position )  % ZX-data 
  else
     pos         = levs; 
  end

  pos(pos==0) = 1;                       % safety to avoid corner pt. problems

%--------------------
% Do the interpolation
%---------------------

  z1   =  zfield(pos+1);
  z2   =  zfield(pos);
  v1   =  var(pos+1);
  v2   =  var(pos);

  tmp  = (z1 - z2) + 1E-35;
  vonz = mask.*( ( (v1-v2)*zval+v2.*z1-v1.*z2 )./tmp);

  return
