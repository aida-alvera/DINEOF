function dataout = cloud_test(datain,mask)

M = size(datain,1);
N = size(datain,2);
Z = size(datain,3);

i = 2:M-1;
j=2:N-1;
   
mask3=repmat(mask,[1 1 Z]);

m = isnan(datain) & mask3==1;
mo = zeros(M,N,Z) == 1;

mo(i,j,:) = m(i+1,j,:) | m(i,j+1,:) | m(i-1,j,:) | m(i,j-1,:) | m(i+1,j+1,:) | m(i-1,j+1,:) | m(i-1,j-1,:) | m(i+1,j-1,:);


dataout=double(mo);
      