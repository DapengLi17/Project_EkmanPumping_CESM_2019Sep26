%% readme 
% constants used in codes
% TAUX is 3600x1800*1440, you can not read such a big matrix in 


%% === no-change constants ===
% infile1 = '../raw_data/CESMncFilesGlobalFromAgartha_ZhangQiuying/CESM_regrid_MonthlyOutput_TAUXTAUYSSH_87-03_test_DL_2020Feb05.nc';
% 
% lon_1d_raw = ncread(infile1,'lon');
% % lon_1d_raw = 0:0.1:359.9;
% lat_1d_raw = ncread(infile1,'lat');
% % lat_1d_raw = -89.95:0.1:89.95;
lon_1d_raw = [0.05   : 0.1 : 360];
lat_1d_raw = [-89.95 : 0.1 : 90]';
g=9.8;
rho_w = 1020;
R = 6371000; % earthRadius('meters')   % Returns 6371000
addpath(genpath('Func4EkmanProject/'))
% ============================


%% === change constants based on your needs ===  
% the original lat is -90 to 90, and original lon is 0 to 360
% lat_limits must be within [-70 70] to discard polar regions. 

% --- Gulf Stream ---
% lat_limits = [-90 90]; % chosen lat ranges for comparison purpose
lat_limits = [25 50]; % chosen lat ranges
lon_limits = [275 325];
step_lat = 1;
step_lon = 1;
% vertical index, 9 elements totally, only use 50m depth for global pic
% same index for temp and wvel
start_z = 1;
count_z = 1;
stride_z = 1; 
% z_w_top = 0, 1000, 2000, 3000, 4000, 5000, 7000, 9000, 11000, [cm]
% z_t   = 500, 1500, 2500, 3500, 4500, 5500, 7500, 9500, 11500, [cm]
% ------------------
% ========================================


%% === compute index based on constants specified above ===
[~,indx4LatMin] = min(abs(lat_1d_raw-lat_limits(1)));% lat_1d_raw(indx4lat_limit1)
[~,indx4LatMax] = min(abs(lat_1d_raw-lat_limits(2)));% lat_1d_raw(indx4lat_limit2)
lat_1d = double(lat_1d_raw(indx4LatMin:step_lat:indx4LatMax));
count_lat = length(lat_1d);

[~,indx4LonMin] = min(abs(lon_1d_raw-lon_limits(1)));% lon_1d_raw(indx4lon_limit1)
[~,indx4LonMax] = min(abs(lon_1d_raw-lon_limits(2)));% lon_1d_raw(indx4lon_limit2)
lon_1d = double(lon_1d_raw(indx4LonMin:step_lon:indx4LonMax));
count_lon = length(lon_1d);

% --- index used for ncread.m ---
%  WVEL   
%            Size:       3600x1800x9x31
%            Dimensions: X,Y,z_w_top,time
%  TEMP
%            Size:       3600x1800x9x31
%            Dimensions: X,Y,z_t,time
start_4dvar  = [indx4LonMin indx4LatMin start_z   1];
count_4dvar  = [count_lon   count_lat   count_z   Inf];
stride_4dvar = [step_lon    step_lat    stride_z  1];
% ------------------------------------------

% %% === find time indx for summer and winter season ===
% % time stamp in CESM
% jultime_raw = [datenum([0087,01,01]):datenum([0087,12,31]) ...
%            datenum([0088,01,01]):datenum([0088,12,31]) ...
%            datenum([0089,01,01]):datenum([0089,12,31]) ...
%            datenum([0090,01,01]):datenum([0090,12,31])];
% jultime_raw(jultime_raw==datenum([0088,02,29]))=[]; % delete 0088 (leap year) Feb 29
% jultime_raw_vec=datevec(jultime_raw);
%   IndxMon = FindMonthlyTimeIndxFunc(jultime_raw_vec);
%   
% % Sea: Season
%   IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
%       IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
%   IndxSea{2} = [IndxMon{4};IndxMon{5};IndxMon{6}; ...
%       IndxMon{7};IndxMon{8};IndxMon{9}];    % Summer season (Jing et al. 2019)
% % ====================================================
