function E = unpack(U,mask);

E = double(mask);
E(find(mask==1)) = U;
E(find(mask==0)) = NaN;