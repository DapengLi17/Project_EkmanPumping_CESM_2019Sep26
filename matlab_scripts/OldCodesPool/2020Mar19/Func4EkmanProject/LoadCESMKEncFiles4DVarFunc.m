%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to load CESM Kuroshio Extension (KE) nc files.
%
% update history:
% v1.0 DL 2019Aug26
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   InFile       - inputfile name (e.g. InFile='kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_TAUX_00870101-00871231.nc')        
%   ParamStr     - parameter string
%                  ParamStr = {'time','ULAT','ULONG','TAUX'} or 
%                  ParamStr = {'time','TLAT','TLONG','SSH'}
% 
% OUTPUT:
%   Time,Lat,Lon,Dat         - structure array, including lat, lon, jultime, temperature (T), 
%                  salinity (S), pressure (P). 
%
% EXTRA NOTES:
% N/A  
% ====================================================================

function [Time,Lat,Lon,Dat,Z] = LoadCESMKEncFiles4DVarFunc(InDir,ParamStr,IndxZ)

%% === load data ===
 file_pattern=fullfile(InDir,['*',ParamStr{1},'*.nc']);
 ncfiles=dir(file_pattern);
 
 ncfile1 = fullfile(InDir,ncfiles(1).name); % the first nc file 
 
 Lon  = ncread(ncfile1,ParamStr{3});
 Lat  = ncread(ncfile1,ParamStr{4});
 Z_raw  = ncread(ncfile1,ParamStr{5});
 Z = Z_raw(IndxZ);
 
  for i = 1 : length(ncfiles)
    clear indx baseFileName inputfile Dat_raw
    baseFileName = ncfiles(i).name;
    inputfile = fullfile(InDir, baseFileName);
    fprintf(1, 'Now reading %s\n', inputfile); 
    indx = [365.*i-364:365.*i];
    Time(indx) = ncread(inputfile,ParamStr{2});
    Dat_raw = ncread(inputfile,ParamStr{1});
    Dat(:,:,:,indx)  = Dat_raw(:,:,IndxZ,:);
  end
 
% ==================

end