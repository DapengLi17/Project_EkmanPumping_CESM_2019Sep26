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

function [Time,Lat,Lon,Dat] = LoadCESMKEncFiles3DVarFunc(InDir,ParamStr)

%% === load data ===
 file_pattern=fullfile(InDir,['*',ParamStr{1},'*.nc']);
 ncfiles=dir(file_pattern);
 
 ncfile1 = fullfile(InDir,ncfiles(1).name); % the first nc file 
 
 Lat  = ncread(ncfile1,ParamStr{4});
 Lon  = ncread(ncfile1,ParamStr{3});
 
  for i = 1 : length(ncfiles)
    clear indx baseFileName inputfile 
    baseFileName = ncfiles(i).name;
    inputfile = fullfile(InDir, baseFileName);
    fprintf(1, 'Now reading %s\n', inputfile); 
    indx = [365.*i-364:365.*i];
    Time(indx) = ncread(inputfile,ParamStr{2});
    Dat(:,:,indx)  = ncread(inputfile,ParamStr{1});
  end
  
% ==================

end
