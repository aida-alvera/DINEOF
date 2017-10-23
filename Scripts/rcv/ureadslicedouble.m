% -----------------------------------------------------
% ------------ GHER file : read function --------------
% ------------ for MATLAB routines --------------------
% ------------ M. Rixen 2000 --------------------------


function tab = ureadslice(fidslice,imax,jmax,kmax,valex,nbmots,nl,ir,pastread)

%ide=1;

c4 = zeros(imax,1);

% need distinction : if imax > nbmots

% load all full records including leading and tailing 4 byte integer

    %data = reshape(fread(fid,nl*(nbmots+2),'single'),nbmots+2,nl);

% remove leading and tailing 4 byte integer

    %c4(1:nbmots*nl) = reshape(data(1:nbmots,:),nbmots*nl,1);

% read last record

    %c4(nbmots*nl+1:end)=fread(fid,ir,'single');

%stop=stopnow

if (nbmots-pastread) > imax
    c4 = fread(fidslice,imax,'double');
    pastread=pastread+imax;
else
    leftinline=nbmots-pastread;
    c4(1:leftinline) = fread(fidslice,leftinline,'double');
    empty=fread(fid,2,'double');
    stilltoread=imax-leftinline;
    c4(leftinline+1:imax) = fread(fidslice,stilltoread,'double');
    pastread=stilltoread;
end


%%%%%%%%%%
  c4(find(c4(:)==valex)) = NaN;
  %tab = reshape(c4,imax,jmax,kmax);
  tab = reshape(c4,imax,1,1);
  
  flag=1;
  
  %fclose(fid);  %modified gzfclose