function X = dineof_pack(Xm,mask);

n = size(Xm,3);
m = sum(mask(:));

X = zeros(m,n); 
  
for k=1:n
    E = Xm(:,:,k)';
    X(:,k) = E(find(mask'==1));
end;


