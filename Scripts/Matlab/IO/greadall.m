function [v,list] = greadall(pattern);

if iscell(pattern)
  list = pattern;
else
  list = filelist(pattern);
end

v = [];
sz = [];

for i=1:length(list)
  disp(['loading ' list{i}]);
  tmp = gread(list{i});
  
  if isempty(v)
    sz = size(tmp);
    v = zeros([prod(sz) length(list)]);
  end
  
  v(:,i) = tmp(:);
end

v = reshape(v,[sz length(list)]);
    
  
