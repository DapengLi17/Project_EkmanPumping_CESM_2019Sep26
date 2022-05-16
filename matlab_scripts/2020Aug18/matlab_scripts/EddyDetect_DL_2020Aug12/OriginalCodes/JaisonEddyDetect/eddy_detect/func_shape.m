function   [match ndim]  =  func_shape ( var, chkdim )
%
% Description : Check the array/matrix shape of the given variable, against 
%               the given number of dimensions. Returns 1 if success, else 0.
%
% Input : var    - (any size/shape) check the shape of this variable
%         chkdim - (scalar)         number of expected dimensions
% 
% Output : match - (scalar) 1 if var is of chkdim size/shape, else 0
%          ndim  - (scalar) number of valid (length > 1) dimensions for var, 
%
%                            0     scalar
%                            1     vector
%                            2     2D matrix
%                            3     3D matrix
%                            4     4D matrix   ......... and so on
%
% Example : a = ones(3,3,2) ;
%           func_shape(a,2)          % --> returns 0, since "a" is of shape/size 3
%           [b c] = func_shape(a,3) ;% --> returns 1, success
% 
%
% Notes  : - Singleton dimensions are not counted to determine match or ndim.
%            Hence, even if match=1, always "squeeze" variables before using 
%            them.
%
% Written By : UCLA ROMS Team (jaison@atmos.ucla.edu)
% Written On : June/11/2008
% Copyright @ UCLA ROMS Team
% Tool       : Eddy Tracker 
% Version    : 1.0
%
%---------------------------------------------------------------------------

     ndim  = sum ( size(var) > 1 ) ;
     match = ndim == chkdim        ;

     return
