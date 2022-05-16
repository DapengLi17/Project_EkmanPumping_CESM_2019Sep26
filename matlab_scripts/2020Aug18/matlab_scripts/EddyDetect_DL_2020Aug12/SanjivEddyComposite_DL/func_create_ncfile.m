function func_create_ncfile(filein)

  nccreate(filein,'fid',   'Dimensions',{'eddynum',Inf},'Datatype','int64','Format','netcdf4')
  nccreate(filein,'tid',   'Dimensions',{'eddynum',Inf})
  nccreate(filein,'x',     'Dimensions',{'eddynum',Inf})
  nccreate(filein,'y',     'Dimensions',{'eddynum',Inf})
  nccreate(filein,'i',     'Dimensions',{'eddynum',Inf})
  nccreate(filein,'j',     'Dimensions',{'eddynum',Inf})
  nccreate(filein,'rkm',   'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'vort',  'Dimensions',{'eddynum',Inf},'Datatype','single')
  nccreate(filein,'slaave','Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'sstave','Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'sssave','Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'sbl'   ,'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'utau'  ,'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'vtau'  ,'Dimensions',{'eddynum',Inf},'Datatype','single')
%DL  nccreate(filein,'qnet'  ,'Dimensions',{'eddynum',Inf},'Datatype','single')

  nccreate(filein,'tauxave'  ,'Dimensions',{'eddynum',Inf},'Datatype','single') %DL
  nccreate(filein,'tauyave'  ,'Dimensions',{'eddynum',Inf},'Datatype','single') %DL 


return 
