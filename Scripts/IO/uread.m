% [flag,data,imax,jmax,kmax,valex,nbmots,format] = uread(file)
% 
% low-level routine for loading Fortran files
%
% Input:
%  file: filename
%
% Output:
%  flag: 1 is successful
%  c4: data
%  imax,jmax,kmax: size of the array
%  valex: exclusion value
%  nbmots: word length
%
% Example:
%
% [flag,c4,imax,jmax,kmax,valex,nbmots] = uread('file.TEM')
%
% M. Rixen 2000
% Alexander Barth, 2007-08-03

% optimisation: 03/04/2003
% Alexander
%
% support of double precision
% Alexander
%
% 2008-07-14
% detecting binary format
% Alexander

function [flag,data,imax,jmax,kmax,valex,nbmots,format] = uread(file)
flag=0;

% all formats to try
formats = {'ieee-be','ieee-le','native'};

for i=1:length(formats)
  fid=gzfopen(file,'r',formats{i});
  if fid==-1
    error(['uread: ' file ' not found.']); 
  end 

  % ignore header

  fread(fid,20,'int32');
  reclen = fread(fid,1,'int32');

  if reclen==24
    break
  else
    fid = gzfclose(fid);  
  end
end

format = formats{i};

%disp(['uread: file is ' formats{i}]);
if fid==0
  error(['uread: format of ' file ' is not recognized.']); 
end 

imax=fread(fid,1,'int32');
jmax=fread(fid,1,'int32');
kmax=fread(fid,1,'int32');
iprec=fread(fid,1,'int32');
nbmots=fread(fid,1,'int32');
valex=fread(fid,1,'float');
fread(fid,2,'int32');

if (iprec==4)
  vtype = 'float';
else
  vtype = 'double';
end

nl=fix((imax*jmax*kmax)/nbmots);
ir=imax*jmax*kmax-nbmots*nl;
ide=1;

if imax<0 | jmax<0 | kmax<0
  nl=0;
  ir=4;
  % disp('Degenerated matrix');
end

% load all full records including leading and tailing 4 byte integer
% for efficiency these integers are read as two floats or one double

if (iprec==4) 
  data = reshape(fread(fid,nl*(nbmots+2),vtype),nbmots+2,nl);
else
  data = reshape(fread(fid,nl*(nbmots+1),vtype),nbmots+1,nl);  
end  

% remove leading and tailing 4 byte integer

data = data(1:nbmots,:);
data = data(:);

% read last record

data = [data; fread(fid,ir,vtype)];

flag=1;
fid = gzfclose(fid);

