function Xm = dineof_unpack(X,mask);

n = size(X,2);
Xm = zeros(size(mask,1),size(mask,2),n);

  
for k=1:n
  Xm(:,:,k) = unpack(X(:,k),mask')';
end


