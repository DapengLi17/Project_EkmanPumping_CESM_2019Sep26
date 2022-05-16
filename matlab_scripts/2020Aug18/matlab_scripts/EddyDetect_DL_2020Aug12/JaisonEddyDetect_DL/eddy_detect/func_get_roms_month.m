function [cday,cmon,cyr,thedate] = get_roms_date_greg(time,tinit,calendar)
 
    % julian
    % gregorian
    % 360_day
    % noleap

     yr0 = tinit(1) ;
     mn0 = tinit(2) ; % not implimented
     dy0 = tinit(3) ; % not implimented


     % ROMS time in seconds
           year=0;
           mois=0;
           day=0;
           month='';
           thedate='';

     ndays = floor(time/(24*60*60)) + 1 ;         % number of days after model initialization
     calendar = lower(calendar) ;
     md_cuml = zeros([13,1]) ;
     nyrs = 0 ;
     if calendar(1:3) == '360'
        md_cuml = md_cuml + [0:30:360]'  ;
        year_len = 360 ;
        nyrs  = int32(fix(ndays/year_len));    % number of years after model init.
     elseif calendar(1:3) == 'jul' | calendar(1:3) == 'gre'
        error('Error:Handle', '\n    FATAL ERROR (get_roms_month) : Only 360_DAY and NOLEAP calendars are allowed. \n')
        md_cuml = [0;31;60;90;121;152;182;213;244;274;305;335;366]  ;
        year_len = 365 ; % current years leap day will be accounted automatically
        nyrs  = int32(fix(ndays/year_len));    % number of years after model init.
     elseif calendar(1:3) == 'nol' 
        md_cuml = [0;31;59;89;120;151;181;212;243;273;304;334;366] ;
        year_len = 365 ;
        nyrs  = floor(ndays/year_len);    % number of years after model init.
     end
 
     cyr      = yr0 + nyrs ;                     % current year 
     cyr_days = ndays - (double(nyrs))*year_len; % days in current year
%display (['   ' num2str(cyr_days)  '   ' num2str(cyr)])

     Month = [ 'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';...
                      'Sep';'Oct';'Nov';'Dec'];
      
     if (cyr_days > md_cuml(13)) 
       error('Error:Handle', '\n    FATAL ERROR (get_roms_month) : Mismatch between specified calendar type and ROMS time. \n')
     elseif (cyr_days == 0)
          cyr = cyr -1 ;
          cyr_days = md_cuml(13)-1;
     end

     mon_prev = md_cuml(md_cuml <= (cyr_days-1)) ;
     cmon     = length(mon_prev) ;
     mname    = Month(cmon,:) ;
     cday     = cyr_days - md_cuml(cmon)  ;
     if (cday > 31 || cday < 1)
       error('\n    FATAL ERROR (func_get_roms_month) : Days in a month is %5i...something is wrong. \n', cday)
     end
     if (cmon > 12 || cmon < 1)
       error('\n    FATAL ERROR (func_get_roms_month) : Months in a year is %5i...something is wrong. \n', cmon)
     end
     thedate  = [num2str(cday),'/',mname,'/',num2str(cyr)] ;
   %  display (['        ',thedate])  
