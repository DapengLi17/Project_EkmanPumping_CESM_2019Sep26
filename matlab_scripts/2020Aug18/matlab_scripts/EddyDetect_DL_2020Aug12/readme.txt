%% === readme ===

% descrip: Matlab codes for eddy detection. The original codes are Jaison's eddy detect Matlab codes for ROMS and Sanjiv's Matlab codes for eddy composites (see OriginalCodes). I used Python (see Project_EkmanPumping_2019Sep27/python_scripts) to convert POP ouput nc files to ROMS-like nc files to use Jaison's Matlab codes. I editted Jaison's and Sanjiv's codes for my personal use (JaisonEddyDetect_DL and SanjivEddyComposite_DL). The changes are listed below. 

1. new_tracker_GulfStream_DL_2020Aug12.m is editted based on Sanjiv's new_tracker_kuro.m (see OriginalCodes/SanjivEddyComposite). Changes in new_tracker_GulfStream_DL_2020Aug12.m are marked with %DL. 

2. func_init_grid.m use Jaison's version since it computes dx, dy, f, pm, pn. Sanjiv's version comment them out. func_init_grid.m line#78-79, replace error() with fprintf(), see below:
%DL            error('\n   FATAL (%s) : Failed to retrieve variables from \n \t %s \n',mfilename, gridfile)
             fprintf('\n   WARNING : Missing values are estimated from func_init_grid.m.\n') %DL


3. func_create_ncfile.m: I comment lines below:
%DL  nccreate(filein,'vort',  'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'sstave','Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'sssave','Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'sbl'   ,'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'utau'  ,'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'vtau'  ,'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'qnet'  ,'Dimensions',{'eddynum',Inf},'Datatype','single')
and add: 
  nccreate(filein,'tauxave'  ,'Dimensions',{'eddynum',Inf},'Datatype','single') %DL
  nccreate(filein,'tauyave'  ,'Dimensions',{'eddynum',Inf},'Datatype','single') %DL 

4.  func_get_hslice.m: I comment lines below:
%DL        hslice.sst  = squeeze( ncr{'temp'}(tindx,zsize,:,:) ) .* grid.mask  ; % it is a 2D var
%DL        hslice.sss  = squeeze( ncr{'salt'}(tindx,zsize,:,:) ) .* grid.mask  ; % it is a 2D var
%DL        hslice.utau = u2rho_2d( squeeze( ncr{'sustr'}(tindx,:,:) ) ) .* grid.mask  ; % it is a 2D var
%DL        hslice.vtau = v2rho_2d( squeeze( ncr{'svstr'}(tindx,:,:) ) ) .* grid.mask  ; % it is a 2D var
%DL        hslice.qnet = squeeze( ncr{'shflux'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var

%DL        hslice.sbl  = squeeze( ncr{'Hsbl'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var

and add:
        hslice.taux = squeeze( ncr{'taux'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var %DL
        hslice.tauy = squeeze( ncr{'tauy'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var %DL

5. func_get_diag_sla_geo.m: I comment lines below:
    % SST
%DL     diag.sst = prog.sst ;

    % SSS
%DL     diag.sss = prog.sss; 

    % UTAU, VTAU, QNET
%DL     diag.utau= prog.utau;
%DL     diag.vtau= prog.vtau;
%DL     diag.qnet= prog.qnet;

    % SBL
%DL     diag.sbl = prog.sbl; 

6. func_get_eddies_circ_102010.m: I use Sanjiv's version (his version comments unnecessary codes in the end compared to Jaison's original code)

7.  

8. 

9.




