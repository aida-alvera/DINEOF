function tab = uload(name)

i = find(name == '#');

if isempty(i)
  [flag,c4,imax,jmax,kmax,valex,nbmots] = uread(name);
  
  %ingreadimax=imax
  %ingreadjmax=jmax
  %ingreadkmax=kmax
  
  if (imax < 0)
    imax = abs(imax);
    jmax = abs(jmax);
    kmax = abs(kmax);    
    
    tab = zeros(imax,jmax,kmax);
    [i,j,k] = ndgrid([1:imax],[1:jmax],[1:kmax]);
    tab = c4(1) + c4(2)*i + c4(3)*j + c4(4)*k;
  else
    
  c4(find(c4(:)==valex)) = NaN;
  tab = reshape(c4,imax,jmax,kmax);
%sizereshapedingread=size(tab)
  end
else
  fname = name(1:i-1);
  vname = name(i+1:end);
  
  zipped = 0;

  if length(fname) > 3
    if strcmp(fname(end-2:end),'.gz') 
       zipped = 1;
    end
  end

  if zipped
    tmp = ['/tmp/gread.' num2str(floor(rand*1e7),'%7.7i')];
%    disp(['!cp ' fname '  ' tmp '.gz;   gunzip ' tmp '.gz; ']);
    eval(['!cp ' fname '  ' tmp '.gz;   gunzip ' tmp '.gz; ']);  
    tab = gread([tmp '#' vname]);
%    disp(['!rm  ' tmp ]);
    eval(['!rm  ' tmp ]);    
  else
  
  nf = netcdf(fname,'nowrite');
  nv = nf{vname};
%  nv = nf{vname,1};
  tab = nv(:);
  tab = permute(tab,[ndims(tab):-1:1]);
  tab = tab(:,:,end:-1:1);
  if (~isempty(nv.missing_value(:)))
    tab(find(tab == nv.missing_value(:))) = NaN;
  end
  close(nf);
  end
  
end