% [flag]=uwrite(file,c4,imax,jmax,kmax,valex,nbmots)
% 
% low-level routine for writing data redable by a Fortran program
%
% Input:
%  file: filename
%  c4: data
%  imax,jmax,kmax: size of the array
%  valex: exclusion value
%  nbmots: word length
%
% Output:
%  flag: 1 is successful
%
% Example:
%
% [flag] = uread('file.TEM',a,size(a,1),size(a,2),size(a,3),9999,1000)
%
% M. Rixen 2000

% -----------------------------------------------------
% ------------ GHER file : write function --------------
% ------------ for MATLAB routines --------------------
% ------------ M. Rixen 2000 --------------------------

%
% 2008-07-14
% change format with global variable UWRITE_DEFAULT_FORMAT
% Alexander
%
% 2010-01-06
% NaNs are no longer substituted by valex since this is done
% in gwrite

function [flag]=uwrite(file,c4,imax,jmax,kmax,valex,nbmots)
global UWRITE_DEFAULT_FORMAT
global UWRITE_DEFAULT_PRECISION

if isempty(UWRITE_DEFAULT_FORMAT)
  format = 'ieee-be';
else
  format = UWRITE_DEFAULT_FORMAT;
end

if isempty(UWRITE_DEFAULT_PRECISION)
  iprec=4;
else
  iprec=UWRITE_DEFAULT_PRECISION;
end

if iprec == 4
  float_precision = 'single';
else
  float_precision = 'double';
end

[file] = gread_tilde_expand(file);

flag=0;
file;
dummy=0;
dummyf=0.;
dummy24=24;

fid=fopen(file,'w',format);
if fid==-1
  error(['uwrite: could not write ' file]);
end

if nbmots==-1
  nbmots=imax*jmax*kmax
end

for i=1:10
  fwrite(fid,dummy,'int32');
  fwrite(fid,dummy,'int32');
end
nl=fix((imax*jmax*kmax)/nbmots);
ir=imax*jmax*kmax-nbmots*nl;
dummyval2=4*nbmots;
fwrite(fid,dummy24,'int32');
fwrite(fid,imax,'int32');
fwrite(fid,jmax,'int32');
fwrite(fid,kmax,'int32');
fwrite(fid,iprec,'int32');
fwrite(fid,nbmots,'int32');
fwrite(fid,valex,'single');
fwrite(fid,dummy24,'int32');

if imax<0 | jmax<0 | kmax<0
  nl=0;
  ir=4;
  disp('Degenerated matrix');
end

ide=1;

for kl=1:nl
  fwrite(fid,4*nbmots,'int32');
  fwrite(fid,c4(ide:ide+nbmots-1),float_precision);
  fwrite(fid,4*nbmots,'int32');
  ide=ide+nbmots;
end

fwrite(fid,4*ir,'int32');
fwrite(fid,c4(ide:ide+ir-1),float_precision);
fwrite(fid,4*ir,'int32');

flag=1;
fclose(fid);


