function ncisz = func_write_eddy_fields(fileout, snap, num, cyc, anti, missing, filein, ncisz)
%
% Description : Write computed eddy tracking fields to output NetCDF file, after
%                  replacing NaN's by missing_values.
%
%
%   snap - current snapshot with respect to input ROMS file, to read
%                ocean_time
%   tid  - current time index with respect to output Eddy file
%
%   tlev - 1 if there is no previous snapshot wrt cyc and anti
%          2 if there is a  previous snapshot wrt cyc and anti
%
%     num  - (struct, three scalars, num.c, num.a, num.m) - number of cyclones 
%                (num.c), number of anticyclones (num.a) for current snap
%                and num.m - number of minimum snapshots (netcdf arrays are 
%                different with respect to this var)
%
%  
%   Independent of min_life, but use func_min_life at the end to get properly
%          shaped NetCDF file. (Jul/02/2008)
%   
%  %   June/23/2009 modified to match the new NetCDF capabilities (tracking is
%         done separately).
%---------------------------------------------------------------

    global vnames 

    matmax = 1.e+15 ;

    tid = num.t ;

    ncilen = max (num.a, num.c) ;


    if ( ncisz < 0  )
        error ([mfilename ':argchk'], '\n\t FATAL (%s) : ncisz should be >= 0 (ncisz=%g). \n', mfilename, ncisz )
    end 


    ncr = netcdf ( filein, 'read' ) ;
       time = ncr{'ocean_time'}(snap) ;
    ncr = close ( ncr ) ;


    nvar   = length(vnames) ;

    nc = netcdf (fileout,'write') ;
        
        % Time

        nc{'time'}(tid) = time ;

        % Eddy ID

        if ( ncilen > ncisz ) 
            nc{'xeddynum'}(ncisz+1:ncilen) = [ncisz+1:ncilen] ;
            err = sync (nc) ;
        end

        % a_ --> anticyclone, c_ --> cyclone, _net = netcdf name & _mat = matlab name
        for iv = 1:nvar ;

           % cyclones  (num.c)

           c_net   = ['cyc_', vnames{iv}] ;
           c_mat   = ['cyc.', vnames{iv}] ;

           eval( [c_mat,'(abs(',c_mat,')>matmax) = missing',';' ]) ;
           eval( [c_mat,'(isnan(',c_mat,')) = missing',';' ]) ;
           nc{c_net}(1:num.c,tid)       = eval( c_mat );  

%           nc{c_net}(1:num.c,tid)       = [1;2];  
%           nc{c_net}(:,:) 

%           pause

           % anticyclones

           a_net   = ['anti_', vnames{iv}] ;
           a_mat   = ['anti.', vnames{iv}] ;
           eval( [a_mat,'(abs(',a_mat,')>matmax) = missing',';' ]) ;
           eval( [a_mat,'(isnan(',a_mat,')) = missing',';' ]) ;
           nc{a_net}(1:num.a,tid)       = eval( a_mat ) ; 
        end

    nc = close (nc) ;

    % adjust ncisz for next round

    ncisz = ncilen ;

    return
