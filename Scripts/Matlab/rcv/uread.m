% -----------------------------------------------------
% ------------ GHER file : read function --------------
% ------------ for MATLAB routines --------------------
% ------------ M. Rixen 2000 --------------------------

% optimisation: 03/04/2003
% Alexander


function [flag,c4,imax,jmax,kmax,valex,nbmots] = uread(file)
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
ide=1;

if imax<0 | jmax<0 | kmax<0
 nl=0;
 ir=4;
% disp('Degenerated matrix');
end

% allocation
%nbmots
%nl
%ir
c4 = zeros(nbmots*nl+ir,1);

% load all full records including leading and tailing 4 byte integer

data = reshape(fread(fid,nl*(nbmots+2),'single'),nbmots+2,nl);

% remove leading and tailing 4 byte integer

c4(1:nbmots*nl) = reshape(data(1:nbmots,:),nbmots*nl,1);

% read last record

c4(nbmots*nl+1:end)=fread(fid,ir,'single');

flag=1;
fclose(fid);  %modified gzfclose
