function mu_2 = mu2(dataini,datarec,mask,dx,dy);
%
% function mu_2 = mu2(dataini,datarec,mask,dx,dy);
% calculation of observational error variance mu^2
% mu^2: variance not retained within the EOF expansion.
%
% input parameters:
%     dataini: initial (cloudy) data set
%     datarec: DINEOF-reconstructed data set
%     mask: land-sea mask used in DINEOF
%     dx, dy: zonal and meridional resolution (in km)
%
% output parameter:
%     mu_2: observational error variance
%
% screen output:
%     several suggested values for mu^2_effective:
%    (mu^2 taking into account spatial correlation of observations)
%     mu^2_eff = mu^2 * r;
%     r = L^2/(dx * dy);
%    (L: correlation length of observational error)
%
% For more info see:  
% J.-M. Beckers, A. Barth, and A. Alvera-Azcarate
%       DINEOF reconstruction of clouded images including 
%       error maps - application to the Sea-Surface Temperature 
%       around Corsican Island. Ocean Science - Volume 2, Number 2
%       pp. 183-199.

sea_1 = sum(mask(:));
land_1 = size(mask,1)*size(mask,2)-sea_1;

sea_all = sea_1*size(dataini,3);
land_all = land_1*size(dataini,3);

mask3d=repmat(mask,[1 1 size(dataini,3)]);
dataini(mask3d==0)=NaN;
missing_all = sea_all - sum(~isnan(dataini(:))) ;

present_all = sea_all - missing_all;

x= dataini(~isnan(dataini));
xr = datarec(~isnan(dataini));

mu_2 = sum(x.^2-xr.^2)/present_all; 

disp([' ']);
disp(['Observational error variance, mu^2 = ' num2str(mu_2)]);
disp([' ']);
disp(['Suggested mu2_eff:']);
disp([' ']);
disp(['for L = 25  km, mu2_eff = ' num2str(25^2*mu_2 /(dx*dy))]);
disp(['for L = 50  km, mu2_eff = ' num2str(50^2*mu_2 /(dx*dy))]);
disp(['for L = 100 km, mu2_eff = ' num2str(100^2*mu_2 /(dx*dy))]);
disp([' ']);
