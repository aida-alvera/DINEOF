% tab = gread(name)
% 
% loads a NetCDF variable from a NetCDF file
%
% Input:
%   name: file name and variable name concanenated by a #
%
% Output:
%   tab: the loaded variable
%
% Example:
%
% Tair = gread('atm.nc#Tair')
%
% Alexander Barth, 2007-08-03

function [tab,info] = gread(name)

[name] = gread_tilde_expand(name);

extraction = [];

% search for extraction list

i = find(name == '(');

if ~isempty(i)
  extraction = name(i:end);
  name = name(1:i-1);
end 


i = find(name == '#');

if isempty(i)
  [flag,c4,imax,jmax,kmax,info.valex] = uread(name);
  if (imax < 0)
    imax = abs(imax);
    jmax = abs(jmax);
    kmax = abs(kmax);    
    
    tab = zeros(imax,jmax,kmax);
    [i,j,k] = ndgrid([1:imax],[1:jmax],[1:kmax]);
    tab = c4(1) + c4(2)*i + c4(3)*j + c4(4)*k;
  else
    tab = reshape(c4,imax,jmax,kmax);

    if (~isempty(extraction)) 
      eval(['tab = tab' extraction ';']);
    end 
  end
else
  fname = name(1:i-1);
  vname = name(i+1:end);
  i = find(vname == '(');
  
  [tmp] = decompress(fname);
  zipped = ~strcmp(tmp,fname);

  if zipped
    [tab,info] = gread([tmp '#' vname extraction]);
    delete(tmp);
  else
    
    nf = netcdf(fname,'r');
    % without auto scaling
    nv = nf{vname};
    % with auto scaling
    %  nv = nf{vname,1};
    
    
    if isempty(nv)
      error(['variable ' vname ' not found in netcdf file ' fname]);
    end
    
    
    if (isempty(extraction))            
      tab = nv(:);
    else
      % reverse extraction list

      delim = [ strfind(extraction,'(') strfind(extraction,',')  strfind(extraction,')')];
      n = length(delim)-1;
      rev = '('; 
      
      for i=0:n-2; 
        rev = [rev extraction(delim(n-i)+1:delim(n-i+1)-1) ','];
      end;   
      
      rev = [rev extraction(delim(1)+1:delim(2)-1) ')'];
      
      eval(['tab = nv' rev ';']);
    end  
    
    tab = permute(tab,[ndims(tab):-1:1]);
    
    if ~isempty(nv.missing_value(:))
      info.valex = nv.missing_value(:);
    elseif  ~isempty(fillval(nv))
       info.valex =  fillval(nv);
    else
      info.valex = NaN;
    end
    
    close(nf);
  end  
end

if isnan(info.valex)
  info.excl = 0;
else
  info.excl = sum(tab(:) == info.valex);
  tab(tab == info.valex) = NaN;
end

