function func_write_eddy_stats(filein,eddy_stats,vlist_comp,fid,tid)

 xtemp = ncread(filein,'x');  
 nstart = max(size(xtemp)); 
 
%nend = max(size(eddy_stats.r)) ; 

% To start writing from 1st place for the first write. 
% Even before writing anything nstart returns '1' (not 0) at the start.
% This is because at the beginning even an empty xtemp has size 0 x 1
% as it has an unlimited netcdf dimension (eddynum)
 if (isempty(xtemp)==1)
   nstart=0;
 end

%pause

 ncwrite(filein,'fid',repmat(fid,[1 eddy_stats.n]),nstart+1);
 ncwrite(filein,'tid',repmat(tid,[1 eddy_stats.n]),nstart+1);

 for iv=1:numel(vlist_comp)

    vlist_id = ['eddy_stats.',vlist_comp{iv}];
    ncwrite(filein,vlist_comp{iv},eval( [vlist_id] ), nstart+1) 

 end 


%pause






return  
