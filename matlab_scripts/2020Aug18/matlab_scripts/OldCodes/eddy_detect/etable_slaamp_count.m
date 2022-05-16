%
%  Description : Bin eddy count according to size/radius.
%  
%   Sep/10/2009
%--------------------------------------------------------------

% USER INPUT

    need_table  = true;
    table       = 'etable_slaamp_count_inf.txt';
    table_title = 'Eddy count binned wrt SLA Amplitude (cm)';

    eddy_file = 'store_eff/etrack_usw51_sla.nc';
    min_snap  = 15 ;
    dt        = 2  ;% days 

    blo       = 1  ;
    bhi       = 15 ;
    bsz       = 1  ;

% INIT
 
    fprintf ('\t NOTE : Reading variables ......\n')  
    cc =  func_init_eddy ('cyc', eddy_file, {'cont','slamax'}) ;
    aa =  func_init_eddy ('anti', eddy_file, {'cont','slamax'}) ;

    % define variables here

    cvar = abs(cc.cont-cc.slamax) ;
    avar = abs(aa.cont-aa.slamax) ;


%=======================================================================
    fprintf ('\t NOTE : Finding IDs of eddies with min life ......\n')  
    cid = func_get_eid_life (cvar, min_snap) ;
    aid = func_get_eid_life (avar, min_snap) ;    


    fprintf ('\t NOTE : Removing NaNs ......\n')  
    camp = cvar(:,cid) ;
    camp = camp(~isnan(camp)) ;

    aamp = avar(:,aid) ;
    aamp = aamp(~isnan(aamp)) ;

%% Binning part
  
    fprintf ('\t NOTE : Cumulative Binning ......\n')  
    blo     = floor(blo/bsz) * bsz ;
    bhi     = ceil(bhi/bsz) * bsz  ;
    nbin    = 1 + (bhi-blo)/bsz;
    
    [acount bins] = func_bin_csum_inf(aamp,nbin,blo,bhi) ;
    [ccount bins] =  func_bin_csum_inf(camp,nbin,blo,bhi) ;
    
% Figures

    fprintf ('\t NOTE : Plotting Figures ......\n')  
    tot = (acount+ccount)/size(cvar,1) ;

    figure(1) ;   plot(bins,tot,'-r')
    figure(2) ;   warning off ; plot(bins,acount./ccount,'-r') ; 
         axis([blo-1 bhi+1 0 1.5]) ; hold on 
         plot ([blo-1 bhi+1],[1 1],'k') ;  
    warning on ;

% output table  

   if ( need_table ) 
     fprintf ('\t NOTE : Writing output Table ......\n')  
     fid = fopen(table, 'w');
     pdir = pwd ;
     mlife = sprintf('%5i',(min_snap-1)*dt) ;
     fprintf(fid, ' TABLE     :  %s  : Min Life = %s days\n',table_title, mlife);
     fprintf(fid, ' Dir       : %s\n',pwd) ;
     fprintf(fid, ' Eddy File : %s\n',eddy_file) ;
     fprintf(fid, ' Note      : Binning is cumulative\n') ;
     fprintf(fid, ' +++++++++++++++++++\n') ;
     fprintf(fid, ' Bin    NCyc   NAnti\n') ;
     fprintf(fid, ' +++++++++++++++++++\n') ;
     tmp = [bins ;ccount ; acount] ;
     fprintf(fid, ' %5.2f  %7i  %7i\n', tmp);
     fclose(fid);     
   end

   eval (['type ' table]) ;
