% flag = gwrite(name,tab)
% 
% save a file in GHER format or a NetCDF variable in NetCDF file
%
% Input:
%   name: name of the GHER file or NetCDF file name and 
%         variable name concanenated by a #
%   tab: the variable to save
%
%
% Output:
%   flag: 1 on success; 0 otherwise
%
% Example:
%
% For GHER format:
% gwrite('data.TEM',temp)
%
% Alexander Barth, 2008-03-19

function flag = gwrite(name,tab)

[name] = gread_tilde_expand(name);

valex = -9999;
tab(isnan(tab)) = valex;  

i = find(name == '#');

if isempty(i)
  [imax] = size(tab,1);
  [jmax] = size(tab,2);
  [kmax] = size(tab,3);
  nbmots =  10*1024*1024; % 10 MB record length

  [flag]=uwrite(name,tab(:),imax,jmax,kmax,valex,nbmots);
else
  fname = name(1:i-1);
  vname = name(i+1:end);

  if (fexist(fname))
    nf = netcdf(fname,'write');      
  else
    nf = netcdf(fname,'clobber');      
  end

  tab = permute(tab,[ndims(tab):-1:1]);

  % re-use dimensions that are already defined in the netCDF file
  
  for i=1:myndims(tab);
    dim{i} = gendim(nf,size(tab,i));
  end
  
  nf{vname} = ncfloat(dim{:});
  
  nf{vname}(:) = tab;
  nf{vname}.missing_value = ncfloat(valex);
  
  close(nf);
end

function ncdim = gendim(nf,len)
  d = dim(nf);
  
  for i=1:length(d)
    if length(d{i})==len
      ncdim = name(d{i});
      return
    end
  end
  ncdim = ['dim' num2str(length(nf)+1,'%3.3i')];
  nf(ncdim) = len;
return

function d = myndims(a)
d = ndims(a);

if (d==2 & size(a,2) == 1)
  d=1;
end
