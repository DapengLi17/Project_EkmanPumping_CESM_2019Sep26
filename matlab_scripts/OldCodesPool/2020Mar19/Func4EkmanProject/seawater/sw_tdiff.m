function tdiff = sw_tdiff(S,T,P)

% SW_TDIFF  thermal diffusivity 
%===========================================================================
% SW_TDIFF  $Revision: 0.0 $  $Date: 1998/01/19 $
%           Copyright (C) Ayal Anis 1998. 
%
% USAGE:  tdiff = sw_tdiff(S,T,P) 
%
% DESCRIPTION:
%    Calculates thermal diffusivity of sea-water. 
%
% INPUT:  (all must have same dimensions)
%   S  = salinity    [psu      (PSS-78) ]
%   T  = temperature [degree C (IPTS-68)]
%   P  = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   tdiff = thermal diffusivity of sea-water [m^2/s] 
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER:  general purpose
% CALLEE:  sw_cp, sw_dens, sw_tcond

%-------------
% CHECK INPUTS
%-------------
if nargin ~= 3
   error('sw_tdiff.m: Must pass 3 parameters ')
end 

% LET sw_dens.m DO DIMENSION CHECKING

%------
% BEGIN
%------

tdiff = sw_tcond(S,T,P)./sw_cp(S,T,P)./sw_dens(S,T,P);

return      
%=========================================================================

