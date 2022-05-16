%% readme 
% constants used in codes


%% === no-change constants ===
lon_1d_raw = [0.05   : 0.1 : 360];
lat_1d_raw = [-89.95 : 0.1 : 90]';
g=9.8;
rho_w = 1020;
R = 6371000; % earthRadius('meters')   % Returns 6371000
% time stamp in CESM
jultime_raw = [datenum([0087,01,01]):datenum([0087,12,31]) ...
           datenum([0088,01,01]):datenum([0088,12,31]) ...
           datenum([0089,01,01]):datenum([0089,12,31]) ...
           datenum([0090,01,01]):datenum([0090,12,31])];
jultime_raw(jultime_raw==datenum([0088,02,29]))=[]; % delete 0088 (leap year) Feb 29
jultime_raw_vec=datevec(jultime_raw);
% ============================


%% === change constants based on your needs ===  
% the original lat is -90 to 90, and original lon is 0 to 360
% lat_limits must be within [-70 70] to discard polar regions. 

% --- Global ---
 lat_limits = [-70 70]; % chosen lat ranges
 lon_limits = [0 360];
 step_lat = 1;
 step_lon = 1; 
% --------------

% --- Kuroshio Extension (KE) Region ---
% lat_limits = [24 48];   % Jing Zhao's KE region
% lon_limits = [120 200]; % Jing Zhao's KE region
% lat_limits = [27 45];   % Zhang Qiuying's original cut for KE region
% lon_limits = [131 169]; % Zhang Qiuying's original cut for KE region
% step_lat = 1;
% step_lon = 1;
% ---------------------------------

% --- Gulf Stream Region ---
% lat_limits = [24 48]; 
% lon_limits = [-90 -20]; % chosen lon ranges, original lon is 0 to 360
% step_lat = 1;
% step_lon = 1;
% --------------------------

% h_deg = 0.1.*step_lat
indx_z = [1 6 9]; % 9 elements totally
% z_w_top = 0, 1000, 2000, 3000, 4000, 5000, 7000, 9000, 11000,
% ========================================


%% === compute index based on constants specified above ===
[~,indx4LatMin] = min(abs(lat_1d_raw-lat_limits(1)));% lat_1d_raw(indx4lat_limit1)
[~,indx4LonMax] = min(abs(lat_1d_raw-lat_limits(2)));% lat_1d_raw(indx4lat_limit2)
lat_1d = double(lat_1d_raw(indx4LatMin:step_lat:indx4LonMax));
count_lat = length(lat_1d);

[~,indx4LonMin] = min(abs(lon_1d_raw-lon_limits(1)));% lon_1d_raw(indx4lon_limit1)
[~,indx4LonMax] = min(abs(lon_1d_raw-lon_limits(2)));% lon_1d_raw(indx4lon_limit2)
lon_1d = double(lon_1d_raw(indx4LonMin:step_lon:indx4LonMax));
count_lon = length(lon_1d);

% --- index used for ncread.m ---
%  TAUX,TAUY,SSH
%            Size:       3600x1800x31
%            Dimensions: lon,lat,time
start_3dvar  = [indx4LonMin indx4LatMin 1];
count_3dvar  = [count_lon   count_lat   Inf];
stride_3dvar = [step_lon    step_lat    1];

%  WVEL   
%            Size:       3600x1800x9x31
%            Dimensions: X,Y,z_w_top,time
%  TEMP
%            Size:       3600x1800x9x31
%            Dimensions: X,Y,z_t,time
start_4dvar  = [indx4LonMin indx4LatMin indx_z(1)           1];
count_4dvar  = [count_lon   count_lat   length(indx_z)      Inf];
stride_4dvar = [step_lon    step_lat    indx_z(2)-indx_z(1) 1];
% ------------------------------------------


%% === convert lat and lon to x and y ===
[x_2d,y_2d] = ConvertLonLat2XY4UnevenGridsFunc(lon_1d,lat_1d);
% =======================================


%% === compute f ===
f = sw_f(lat_1d); % f = sw_f(lat)
f_rp = repmat(sw_f(lat_1d),1,length(lon_1d));
f_min = sw_f(10); % f = sw_f(lat), min f to discard Ekman Wvel near equator
% ==================


%% === find time indx for summer and winter season ===
  IndxMon = FindMonthlyTimeIndxFunc(jultime_raw_vec);
  
% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
      IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
  IndxSea{2} = [IndxMon{4};IndxMon{5};IndxMon{6}; ...
      IndxMon{7};IndxMon{8};IndxMon{9}];    % Summer season (Jing et al. 2019)
% ====================================================
