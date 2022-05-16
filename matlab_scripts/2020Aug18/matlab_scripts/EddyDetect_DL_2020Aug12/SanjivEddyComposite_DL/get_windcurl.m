function [taucurl]=get_windcurl(taux,tauy,pm,pn,mask)

% The order needs to be (nz,ny,nx) because the older ROMS routines were written assuming variables are read in that way.
pm=pm'; pn=pn'; mask=mask';


if (length(size(taux))==3)
  dimorder=[3 2 1];
elseif (length(size(taux))==2)
  dimorder=[2 1];
end 

%taux=permute(taux,dimorder); tauy=permute(tauy,dimorder); 

mask(2:end-1,2:end-1)=mask(1:end-2,1:end-2).*...
                      mask(1:end-2,2:end-1).*...
                      mask(1:end-2,3:end).*...
                      mask(2:end-1,1:end-2).*...
                      mask(2:end-1,2:end-1).*...
                      mask(2:end-1,3:end).*...
                      mask(3:end,1:end-2).*...
                      mask(3:end,2:end-1).*...
                      mask(3:end,3:end);


if (length(size(taux))==3)

  parfor it=1:size(taux,1);
     taucurl(it,:,:)=mask.*psi2rho(vorticity(squeeze(taux(it,:,:)),squeeze(tauy(it,:,:)),pm,pn));
  end
 
elseif (length(size(taux))==2)

     taucurl(:,:)=mask.*squeeze( psi2rho(vorticity(squeeze(taux(:,:)),squeeze(tauy(:,:)),pm,pn)) );

else

  'Error! The winds need to have at least 2 dimensions!' 

end 

taucurl=permute(taucurl,dimorder); 





