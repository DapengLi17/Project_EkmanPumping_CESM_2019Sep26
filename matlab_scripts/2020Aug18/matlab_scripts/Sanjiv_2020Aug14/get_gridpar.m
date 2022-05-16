function [grid_struc]=get_gridpar(gridfile)


grid_struc.ocean_time = squeeze(ncread(gridfile,'ocean_time')); 
grid_struc.Cs_r    = ncread(gridfile,'Cs_r');
grid_struc.theta_s = ncread(gridfile,'theta_s');
grid_struc.theta_b = ncread(gridfile,'theta_b');
grid_struc.hc      = ncread(gridfile,'hc');
grid_struc.vtrans  = ncread(gridfile,'Vtransform');
grid_struc.vstret  = ncread(gridfile,'Vstretching');

grid_struc.s_rho   = ncread(gridfile,'s_rho');
grid_struc.s_w     = ncread(gridfile,'s_w');
grid_struc.nz_rho  = getfield(ncinfo(gridfile,'s_rho'),'Size');
grid_struc.nz_w    = getfield(ncinfo(gridfile,'s_w'),'Size');

grid_struc.bathy   = ncread(gridfile,'h'); 
grid_struc.pm      = ncread(gridfile,'pm');
grid_struc.pn      = ncread(gridfile,'pn');
grid_struc.lonr    = ncread(gridfile,'lon_rho');
grid_struc.latr    = ncread(gridfile,'lat_rho');
grid_struc.maskr   = squeeze(ncread(gridfile,'mask_rho'));

















end 
