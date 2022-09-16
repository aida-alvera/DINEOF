function dineof_cvp(fname,maskfname,outdir,nbclean)
%
%function dineof_cvp(fname,maskfname,outdir,nbclean)
% Extracts clouds from file 'fname' and adds them to the 'nbclean'
%   cleanest images from the data set
%
% Input variables:
%       fname (string): file name in the disk (gher format)
%       maskfname (string): mask file name in teh disk (gher format)
%       outdir (string): directory where new file will be written
%       nbclean (integer): number of cleanest images to be covered with clouds
%  
% File with added clouds will be written in a file with the same name
%    as the initial file (and at the same location in the disk)
%    but with the extension '.clouds' added. An additional file, 
%    clouds.index is also written, which contains the indexes of 
%    the added clouds (can be used as cross-validation points in 
%    DINEOF)
%
% Reference:
% A. Alvera-Azcarate, A. Barth, M. Rixen, and J. M. Beckers. 
%    Reconstruction of incomplete oceanographic data sets using 
%    Empirical Orthogonal Functions. Application to the Adriatic 
%    Sea Surface Temperature. Ocean Modelling, 9:325-346, 2005.
%


SST = gread(fname);
mask = gread(maskfname);

for k=1:size(SST,3)
  tmp = SST(:,:,k);
  tmp(mask == 0) = NaN;
  SST(:,:,k) = tmp;
end

%rand('state',0);

nbland = sum( mask(:) == 0 );
m = sum( mask(:) == 1 );

cloudcov = (sum(sum(isnan(SST),2),1) - nbland)/m;
cloudcov = cloudcov(:);

%clean = find(cloudcov < 0.2);

[c,clean] = sort(cloudcov);
clean = clean(1:nbclean);

N = length(cloudcov);

index = floor(N * rand(nbclean,1))+1;

while (any(index == clean))
  index = floor(N * rand(nbclean,1))+1;
end

SST2 = SST;
SST2(:,:,clean) =  SST(:,:,clean) ./ ~isnan(SST(:,:,index));
SST2(isinf(SST2)) = NaN;


%mask3d = mask(:,:,ones(size(SST,3),1));

imax = size(SST,1);
jmax = size(SST,2);

m=0;

for i=1:imax    
  for j=1:jmax
    if (mask(i,j) == 1)
      m = m+1;
      mindex(i,j) = m;
      iindex(m) = i;
      jindex(m) = j;
    end 
  end 
end 



%indexex = find(isnan(SST2) & ~isnan(SST)  & mask3d == 1);
indexex = find(isnan(SST2) & ~isnan(SST));
[iex,jex,kex] = ind2sub(size(SST2),indexex);

clouds_indexes = zeros(length(iex),1);

for l=1:length(iex)
  clouds_indexes(l,1) = mindex(iex(l),jex(l));
  clouds_indexes(l,2) = kex(l);
end  

cloudcov2 = (sum(sum(isnan(SST2),2),1) - nbland)/m;
cloudcov2 = cloudcov2(:);

nbgood = sum(~isnan(SST(:)));
nbgood2 = sum(~isnan(SST2(:)));

[num2str(100*(nbgood-nbgood2)/nbgood) ' % of cloud cover added']

clear SST* clear mask* jex jindex iindex indexex iex kex mindex tmp

%gwrite([fname '.clouds'],SST2);
gwrite([outdir '/clouds.index'],clouds_indexes);

