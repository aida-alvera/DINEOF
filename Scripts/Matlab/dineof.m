% filled = dineof(data,...)
% run DINEOF from Matlab/Octave. data is a 3-d matrix with the dimensions 
% longitude, latitude and time. 
%
% INSTALLATION: place the binary dineof (or dineof.exe under PC) in the parent
% directory (relative to this script dineof.m) or set the option 'exec' 
% described below.
%
% Optional input arguments are:
%
% mask: land-sea mask, if it is not provided it is determined by the script 
%   dineof_mask
% frac: minium fraction of data in time to that a pixel is considered sea
%   (default 0.05, i.e. 5 %)
% exec: executable of DINEOF (default: the full path of the parent directory
%   followed by "dineof" or "dineof.exe" (on PC))
% dir: working directory (default a temporary directory)
% 
% For all other optional parameters, please refer to the DINEOF documentation.
% nev
% neini
% tol
% toliter
% rec
% eof
% norm
% seed
% alpha
% numit
% time
% 
% Note this scripts unset the variable LD_LIBRARY_PATH
% te prevent errors likes these in matlab:
% 
% /opt/matlab-R2013a/sys/os/glnxa64/libgfortran.so.3: version `GFORTRAN_1.4' not found (required by /home/abarth/src/dineof/dineof)

function [filled,info] = dineof(data,varargin)

global DINEOF_EXEC

% default values

nev = 10;
neini = 1;
ncv = nev + 10;
tol = 1.0e-8;
nitemax = 300;
toliter = 1.0e-3;
rec = 1;
eof = 1;
norm = 0;
seed = 243435;
alpha = 0.;
numit = 3;
frac = 0.05;
dir = tempname;
time = [];
mask = [];

exec = DINEOF_EXEC;
if isempty(exec)  
  exec = getenv('DINEOF_EXEC');
  if isempty(exec)  
    % assume the dineof binary is in the parent directory
    if ispc
      exec = fullfile(fileparts(mfilename('fullpath')),'..','dineof.exe');
    else
      exec = fullfile(fileparts(mfilename('fullpath')),'..','dineof');
    end
  end 
end

prop = varargin;

for i=1:2:length(prop)
  if strcmp(prop{i},'nev')
    nev = prop{i+1};
  elseif strcmp(prop{i},'neini')
    neini = prop{i+1};
  elseif strcmp(prop{i},'tol')
    tol = prop{i+1};
  elseif strcmp(prop{i},'toliter')
    toliter = prop{i+1};
  elseif strcmp(prop{i},'rec')
    rec = prop{i+1};
  elseif strcmp(prop{i},'eof')
    eof = prop{i+1};
  elseif strcmp(prop{i},'norm')
    norm = prop{i+1};
  elseif strcmp(prop{i},'seed')
    seed = prop{i+1};
  elseif strcmp(prop{i},'alpha')
    alpha = prop{i+1};
  elseif strcmp(prop{i},'numit')
    numit = prop{i+1};
  elseif strcmp(prop{i},'dir')
    dir = prop{i+1};
  elseif strcmp(prop{i},'exec')
    exec = prop{i+1};
  elseif strcmp(prop{i},'time')
    time = prop{i+1};
  elseif strcmp(prop{i},'mask')
    mask = prop{i+1};
  elseif strcmp(prop{i},'frac')
    frac = prop{i+1};
  else
    error(['unkown argument ' prop{i}]);
  end 
end

fprintf('Using dineof binary: %s\n',exec);

if isempty(mask)
  mask = dineof_mask(data,frac);
end

d = [dir filesep];
mkdir(d);

gwrite([d 'data'],data);
gwrite([d 'mask'],mask);

if alpha > 0
  gwrite([d 'time'],time);
end



initfile = [d 'dineof.init'];
fid = fopen(initfile,'w');

fprintf(fid,'# inpput File for dineof\n');
fprintf(fid,'# created by the script %s\n',mfilename);
fprintf(fid,'data    = [''%s'']\n',[d 'data']);
fprintf(fid,'mask    = [''%s'']\n',[d 'mask']);
fprintf(fid,'nev     = %d\n',nev);
fprintf(fid,'neini   = %d\n',neini);
fprintf(fid,'ncv     = %d\n',ncv);
fprintf(fid,'tol     = %g\n',tol);
fprintf(fid,'nitemax = %d\n',nitemax);
fprintf(fid,'toliter = %g\n',toliter);
fprintf(fid,'rec     = %d\n',rec);
fprintf(fid,'eof     = %d\n',eof);
fprintf(fid,'norm    = %d\n',norm);
fprintf(fid,'Output  = ''%s''\n',[d]);
fprintf(fid,'results = [''%s'']\n',[d 'output']);
fprintf(fid,'seed    = %d\n',seed);

fprintf(fid,'EOF.U = [''eof.nc#U''] \n');
fprintf(fid,'EOF.V = ''eof.nc#V''\n');
% must be there but it is not used
fprintf(fid,'EOF.Sigma = ''eof.nc#Sigma''\n');


fprintf(fid,'alpha   = %g\n',alpha);
fprintf(fid,'numit   = %d\n',numit);

if alpha > 0
  fprintf(fid,'time    = ''%s''\n',[d 'time']);
end

fclose(fid);

cmd = sprintf('"%s" "%s"',exec,initfile);
disp(cmd)

cmd = ['unset LD_LIBRARY_PATH; '  cmd];

[status,result] = system(cmd,'-echo');

if status ~= 0
  disp(result)
  error('dineof failed');
end

filled = gread([d 'output']);

ncfile = fullfile(d,'eof.nc');
info.U = ncread(ncfile,'U');
info.V = ncread(ncfile,'V');
%info.Sigma = ncread(ncfile,'Sigma');
info.Sigma = ncread(fullfile(d,'DINEOF_diagnostics.nc'),'vlsng');
