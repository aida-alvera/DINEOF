function mask = dineof_mask(data,frac);

if ischar(data)
  data = gread(data);
end

count = sum(~isnan(data),3)/size(data,3);
mask = count > frac;

