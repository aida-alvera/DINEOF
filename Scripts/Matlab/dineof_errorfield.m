function [error_field,datarecOI,indi_field]=dineof_errorfield(dataini,datarec,mask,crossval,mu2_eff,vlsng,lftvec,N)
%
%function [error_field,datarecOI]=dineof_errorfield(folder,datafname,resultsfname,maskfname,outdir)
%
% DINEOF error field estimations
% To be used after DINEOF is finished. 
%
% Input variables: 
%     dataini: initial (cloudy) data
%     datarec: field reconstructed by DINEOF 
%     mask: mask (0 for land points, 1 for sea points)
%     crossval: indexes the cross-validation points 
%     mu2_eff: error variance corrected by the spatial correlation
%           of the data (see reference)
%     vlsng: singular values of the data, calculated by DINEOF
%              (file outputEof.vlsng)
%     lftvec: Left EOFs of the data, caculated by DINEOF
%              (file outputEof.lftvec)
%     N: number of EOFs retaind for the DINEOF reconstruction
%              (file neofretained.val)
%     
%
% Output variables:
%     err_field: error field (same units as initial data) for each time step
%           of the data set
%     datarecOI: data set reconstructed using an EOF-based OI (see reference)
%     outliers: flag field. ranging from 0 (perfect) via 3 (suspect data), to 6 and above (very suspect data)
%    
%
%
% Reference:
% J.-M. Beckers, A. Barth, and A. Alvera-Azcarate
%       DINEOF reconstruction of clouded images including 
%       error maps - application to the Sea-Surface Temperature 
%       around Corsican Island. Ocean Science - Volume 2, Number 2
%       pp. 183-199.


mu2_values = mu2_eff;

%Singular values
  Sigmat=vlsng; 
%Left EOFs
  Ut=lftvec;
  clear lftvec vlsng

  mask(find(isnan(mask)))=0;

%Cross-validation points
  cvp_indexes = crossval;
  clear crossval

% number of time instants
   n=size(dataini,3);

% number of sea points (clouded or not)
% equal to size(Ut,1)

  m = sum(mask(:));

%Retained EOFs
  U=Ut(:,1:N); 
  Sigma = Sigmat(1:N);

% scaled EOFs
  L = 1/sqrt(n) * U * diag(Sigma);

% Ut and Vt are strored in C-style 
% transpose the mask and matrix !

  X = zeros(m,n);    
  for k=1:n; 
    X(:,k) = dineof_pack(dataini(:,:,k),mask);   
  end;  
  clear dataini

%Total missing data points
  M = sum(isnan(X(:)));



  meanX =  nanmean(X(:));
  
  stdX = 1;

  % demean 

  X = (X-meanX)/stdX;

  if ~isempty(cvp_indexes)
  % isolate cross-validation points
  icross = size(cvp_indexes,1);
  whos X cvp_indexes
%cross-validation points indexes
  cvp_linindex = sub2ind(size(X),cvp_indexes(:,1),cvp_indexes(:,2));  

%extract cross-validation points from initial matrix 
  Xcvp = X(cvp_linindex) ;
%and take them out of the initial matrix
  X(cvp_linindex) = NaN;
   
  end

%L = 1/sqrt(n) * U * diag(Sigma);

% mask of present data
% 1 (true) means present
% 0 (false) means missing

present = ~isnan(X);

var = sum(L.^2,2);
Xr = dineof_pack(datarec,mask);

for mu2 = mu2_values

  % d are data present, absent data are zero
  d = X;
  d(isnan(d)) = 0;

  %ntime = 1;
  % loop over all time instantes
  err = zeros(m,n);
  err2 = zeros(m,n);
  relerr = zeros(m,n);
  Xr2 = zeros(m,n);
  Xr2p = zeros(m,n);
  sqrtC = zeros(N,N,n);
  
  % sqrt of error covariance  of a particular time instance
  % Es * Es' is error coviarance E
  Es = zeros(m,N);
  err_sum = zeros(n,1);

%  all_ntimes = unique(cvp_indexes(:,2));

% only cross-validation time instants  
%  for i_ntime = 1:length(all_ntimes)
%    ntime = all_ntimes(i_ntime);

  for ntime=1:n
    disp(['Calculating at ntime = ' num2str(ntime)]);
    
    %Lp = L(find(present(:,ntime)),:);

    % or better:
    Lp = L;
    Lp(~present(:,ntime),:) = 0;

    Cp = Lp'*Lp;
    invD =  inv(Cp + mu2 * eye(N,N));
    %C = eye(N,N) -  invD * Cp;
    C = mu2 *  invD;
    
    Xr2(:,ntime) = L * invD * (Lp' * d(:,ntime));
    %  Xr2p(:,ntime) = L * inv(Cp + mu2p * eye(N,N)) * Lp' * d(:,ntime);
 
    if 1
   
    % err2 is epsilon^2 of formula 34
    % this is equivalent and faster

    % old fashioned fortran-style loops

    %for i=1:m
    %  err2(i) = L(i,:) * C * L(i,:)';
    %end


      err2 = sum((L * C) .* L,2);
      err2(err2 < eps) = eps;
      err(:,ntime) = stdX * sqrt(err2);
      relerr(:,ntime) = sqrt(err2./var);
      
      sqrtC(:,:,ntime) =  chol(C)';
      Es = L * sqrtC(:,:,ntime) ;
% JMB
      Esp = Lp * sqrtC(:,:,ntime) ;
           
%DELTASQUARED= mu2 * ones(N,1)- sum(Esp.*Esp,2);
%      outliers_m=zeros(length(Xr),1);
      DELTASQUARED= mu2 * ones(m,1)- sum(Esp.*Esp,2);
      DELTASQUARED=sqrt(DELTASQUARED);
      outliers(:,ntime)=abs((X(:,ntime)-Xr2(:,ntime))./DELTASQUARED(:));      
      outliers_S = nanmedian (outliers(:,ntime));
      outliers_d = 1.4826 * nanmedian(abs(outliers(:,ntime)-outliers_S));
      indi(:,ntime)=abs(outliers(:,ntime)-outliers_S)/outliers_d;
      
 %     outliers_m(abs(outliers(:,ntime)-outliers_S)>M_par*outliers_d)=1;
      
      
% mettre encore un 0 ou NaN aux points sans données
%      outliers(:,ntime)=outliers(:,ntime).*present(:,ntime);
%      outliers(:,ntime)=outliers_d.*present(:,ntime);
       
% end JMB
      err_sum(ntime) = sum(sum(Es.^2,1));
    end
  end

  % error std-dev of spatial average taking the covariance into account
  
  err_sum = sqrt(err_sum)/m;

  % mean error under clouds
  mean_err =  sqrt(mean(err(~present).^2));

  % mean error everywhere
  
  mean_err_global =  sqrt(mean(err(:).^2));
   
  % mean RMS difference between DINEOF and EOF-based OI
  Xr = dineof_pack(datarec,mask);
  Xr = Xr - meanX;
  rms_oi_dineof = sqrt(mean((Xr(:) - Xr2(:)).^2));
  
  % calculate cross validator
  if ~isempty(cvp_indexes)
  val = mean((Xcvp - Xr2(cvp_linindex)).^2);
  val = sqrt(val);
  end
  
  disp([' ']);
  disp(['Mean of initial data, ' num2str(meanX)]);
  disp(['Mean of DINEOF reconstructed data, ' num2str(nanmean(Xr(:)) + meanX)]);
  disp(['Mean of EOF-based OI reconstructed data, ' num2str(nanmean(Xr2(:)) + meanX)]);
  disp([' ']);
  disp(['Mean error of the reconstruction (all points), ' num2str(mean_err_global)]);
  disp(['Mean error of the reconstruction (under clouds), ' num2str(mean_err)]);
  disp([' ']);
  if ~isempty(cvp_indexes)
  disp(['RMS between DINEOF and EOF-based OI at cross-validation points, ' num2str(val)]);
  end
  disp(['RMS between DINEOF and EOF-based OI (all values), ' num2str(rms_oi_dineof)]);
 
  
    error_field = zeros(size(mask,1),size(mask,2),n);
    datarecOI = zeros(size(mask,1),size(mask,2),n);
    

    error_field = dineof_unpack(err,mask);
    datarecOI = dineof_unpack(Xr2,mask)+meanX;
%JMB
    indi_field = zeros(size(mask,1),size(mask,2),n);
    indi_field = dineof_unpack(indi,mask);
%ENDJMB
    
  

end
