function track = func_track_eddy_20080102_rpr  ( etype, min_snap, tmpnc )
%
%  function func_track_eddy_20080102  ( etype, min_snap )
%
%  Description : Modified (2008/01/02) eddy "tracking" algorithm.
% 
%
%  Inputs : 1. etype    - character (scalar) - 'cyc' or 'anti', eddy type
%           2. min_snap - float (scalar) - minimum number of snapshots,  
%                             required for a track, should be >= 2.
%           3. tmpnc    - character (scalar) - name of tmp nc file.
%                           (should be a blank copy of tmpnc, with no values
%                            written).
%
%
%  Output : track : float, struct : ts, te, isz, lsz, n, lind, dkm
%
%  NOTE : 1. Writing to the NetCDF file is not feasible within this function,
%              as the dimension xeddynum is defined based on the maximum 
%              number of tracks among anticyclones and cyclones. Since we 
%              have to keep xeddynum as the unlimited dimension for the tmpnc,
%              we cannot have two different xeddynum dimensions for cyc and 
%              anti.
%
%
%  Jan/02/2008
%
%
%  Jan/16/2009 : Implimented the logic used in func_uniq_eddy for tracking:
%                  set the maximum distance between eddies as rad1+rad2.
%                  With circle fitting, this is a very robust test. 
%
%  aug/20/2009 : fixed a small bug with possible eddy spliting scenarios.
%                   (same eddy being attached to different tracks).
%
%  Aug/25/2009 : Added extra criterion than the maximum distance/travel test.
%                   Amplitude change, Area Change
%------------------------------------------------------------------------------


% Initialize
 
    fprintf ('\n\t  Tracking %s .......\n', upper(etype) ) 

    track.ts    = [] ;
    track.te    = [] ;
    track.lind  = [] ;
    track.isz   = [] ;
    track.tsz   = [] ;
    track.n     = [] ;
    track.dkm   = [] ;

    min_ratio   = 0.25 ; % min value for amp and area ratio
    max_ratio   = 2.50 ; % max value for    "

    matmax = 1.e+15 ;

    nc = netcdf ( tmpnc, 'read' ) ;
        ref      =  nc{[etype '_i']}(:)   ;
        eddy.xkm =  nc{[etype '_xkm']}(:)   ;
        eddy.ykm =  nc{[etype '_ykm']}(:)   ;
        eddy.rkm =  nc{[etype '_rkm']}(:)   ;
        eddy.cont =  nc{[etype '_cont']}(:)   ;
        eddy.slamax =  nc{[etype '_slamax']}(:)   ;
        eddy.area=  nc{[etype '_area']}(:)   ;
    nc = close (nc)                ;


    ref(abs(ref) > matmax)           = NaN   ;
    eddy.xkm(abs(eddy.xkm) > matmax) = NaN   ;
    eddy.ykm(abs(eddy.ykm) > matmax) = NaN   ;
    eddy.rkm(abs(eddy.rkm) > matmax) = NaN   ;
    eddy.area(abs(eddy.area) > matmax) = NaN   ;
    eddy.cont(abs(eddy.cont) > matmax) = NaN   ;
    eddy.slamax(abs(eddy.slamax) > matmax) = NaN   ;
    eddy.amp = abs(eddy.cont-eddy.slamax) ;

    [isize tsize] = size(ref) ;
 

    % remove the padded missing values at the end (this has been done to have a common
    %     X-dimension/axis for the eddy tracking output NetCDF file).    

    clear ref_0 ;
    ref_0 = ref ;
    ref_0(isnan(ref_0)) = 0 ; % set NaN values to zero, for reference variable
    isz = find ( sum(ref_0,2),1,'last' ) ; % Find the I/column-index for the last valid eddy.
    if ( isz < isize )
        ref       = ref(1:isz,:)      ;       % re-size all variables used below
        eddy.xkm  = eddy.xkm(1:isz,:) ;       %           "
        eddy.ykm  = eddy.ykm(1:isz,:) ;       %           "
        eddy.rkm  = eddy.rkm(1:isz,:) ;       %           "
        eddy.amp  = eddy.amp(1:isz,:) ;       %           "
        eddy.area = eddy.area(1:isz,:) ;       %           "
        ref_0       = ref_0(1:isz,:) ;
        isize = isz          ; 
    end

    track.isz = isize ; % need later, for accurate linear-indices
    track.tsz = tsize ; % need later, for accurate linear-indices
    
% Extract tracks - linear indices 

    trID = 0 ; % track ID

    for it = 1:tsize-1 ; % iterate through time
    %for it = 1:1 ; % iterate through time
         % determine the id-start and id-end for valid eddies for this snapshot     
         eids = find(ref_0(:,it),1,'first') ; % this will be 1 always
         eide = find(ref_0(:,it),1,'last') ;
         %fprintf('\t Time = %4i EidS = %2i EidE = %2i \n',it, eids, eide) 

         for ie = eids:eide % iterate through each identified eddies 
              if ( ref_0(ie,it) ~= 0 )   % if it is a valid and non-tracked eddy
                 ts = it ;               % starting time point for this eddy
                 tid = 1 ;               % count number of snapshots
                 lind = [] ;             % init lind variable
                 dist = [] ;             % distance travelled by eddy from previous snap
                 xref = eddy.xkm(ie,it) ; 
                 yref = eddy.ykm(ie,it) ;
                 rref = eddy.rkm(ie,it) ;
                 aref = eddy.amp(ie,it) ;
                 ARref= eddy.area(ie,it) ;  
                 lind(tid) = (it - 1 ) * isize + ie ; % find linear indices
                 dist(tid) = 0 ;                      % dkm=0 at very first snap for an eddy
                % fprintf ('\t  Time = %4i  Eid = %4i   X = %6.2f   Y = %6.2f\n',it,ie,xref,yref)
                 for its = it+1:tsize  % from next snap till end
                    [dist_tmp indx_tmp] = func_NNR_indx_dist ( eddy.xkm(:,its), eddy.ykm(:,its), xref, yref) ; % find the NNR at each snapshot starting from next
                    %if ( dist_tmp <= min( eddy.rkm(indx_tmp,its), rref ) ) % found its next position
  
                    amp_ratio = eddy.amp(indx_tmp,its)/aref ;
                    area_ratio = eddy.area(indx_tmp,its)/ARref ;
                    chk_min = min ([amp_ratio, area_ratio, 1/amp_ratio, 1/area_ratio] ) ;
                    chk_max = max ([amp_ratio, area_ratio, 1/amp_ratio, 1/area_ratio] ) ;
                    if ( dist_tmp <= eddy.rkm(indx_tmp,its)+rref && ( chk_min > min_ratio && chk_max < max_ratio ) ) % found its next position, Jan/16/2009, Aug/25/2009
                       
                       tid = tid + 1 ;  
                       dist(tid) = dist_tmp ;
                       indx(tid) = indx_tmp ;  
                       ref_0(indx(tid),its) = 0 ; % mark this as a tracked eddy
                       lind(tid) = (its - 1 ) * isize + indx(tid) ; % linear indices
                       xref = eddy.xkm(indx(tid),its) ; % update ref location for
                       yref = eddy.ykm(indx(tid),its) ; %    next snapshot
                       rref = eddy.rkm(indx(tid),its) ;
                       aref = eddy.amp(indx(tid),its) ;
                       ARref= eddy.area(indx(tid),its) ;

                       % to avoid possible double counting, set the alreay used eddy to NaN values Aug/20/2009
                       eddy.xkm(indx(tid),its) = NaN ;
                       eddy.ykm(indx(tid),its) = NaN ;
                       eddy.rkm(indx(tid),its) = NaN ;
                    else 
                       break ; 
                    end ;  
                 end
                 if ( tid >= min_snap ) 
                     trID           = trID + 1         ; % increment track ID by 1
                     track.lind(trID,:)  = zeros(1,tsize)   ; % init linear indx var
                     track.dkm(trID,:)   = zeros(1,tsize)   ; % init linear indx var
                     track.ts(trID) = ts               ; % track start time
                     track.te(trID) = ts + tid - 1     ; % track end  time 
                     track.lind(trID,1:tid) = lind(:)       ; % linear indx var
                     track.dkm(trID,1:tid)  = dist(:)       ; % distance travelled in km
                   %  fprintf ('\t  Time = %4i  Eid = %4i   X = %6.2f   Y = %6.2f\n',its,indx(tid),xref,yref)
                 end 
             end
         end
    end


    % do it for last time record, if min_life == 1

    if ( min_snap == 1 )
        it = tsize ;
        eids = find(ref_0(:,it),1,'first') ; % this will be 1 always
        eide = find(ref_0(:,it),1,'last') ;
            
        for ie = eids:eide % iterate through each identified eddies 
            if ( ref_0(ie,it) ~= 0 )   % if it is a valid and non-tracked eddy
               trID = trID + 1 ;
               tid  = 1 ;
               track.lind(trID,:)  = zeros(1,tsize)   ; % init linear indx var
               track.dkm(trID,:)   = zeros(1,tsize)   ; % init linear indx var
               track.ts(trID) = it                    ; % track start time
               track.te(trID) = it                    ; % track end  time 
               track.lind(trID,1:tid) = (it - 1 ) * isize + ie        ; % linear indx var
               track.dkm(trID,1:tid)  = 0                          ; % distance travelled in km

            end
        end
    end

    track.n = trID ; % total number of tracks

% DONE 

  return
