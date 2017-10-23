function u = greadfmt(dir_fmt,range)

u = [];
for i=range(:)'
  fname = sprintf(dir_fmt,i);
  disp(fname);

  tmp = gread(fname);
  
  if isempty(u)
    sz = size(tmp);
    u = zeros([prod(sz) length(range)]);
  end

    
  u(:,i) = tmp(:);
end

u = reshape(u,[sz length(range)]);
