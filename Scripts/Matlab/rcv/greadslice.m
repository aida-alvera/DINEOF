function [flag,imax,jmax,kmax,valex,nbmots,nl,ir,pastread,fid] = greadslice(file)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%here under direct reading of file header  
   
flag=0;
imax=0;
jmax=0;
kmax=0;
valex=0;

fid=fopen(file,'r','ieee-be'); %modified ori: gzfopen and ieee.be
if fid==-1
  error([file ' not found.']); 
end 

% ignore header

fread(fid,21,'int32')

imax=fread(fid,1,'int32')
jmax=fread(fid,1,'int32')
kmax=fread(fid,1,'int32')
iprec=fread(fid,1,'int32')
nbmots=fread(fid,1,'int32')
valex=fread(fid,1,'single')
fread(fid,2,'int32');

nl=fix((imax*jmax*kmax)/nbmots);
ir=imax*jmax*kmax-nbmots*nl;
pastread=0;

flag=1;
