% fid=gzfopen(fname,permisson,machineformat)
% 
% same as fopen except that the file is decompressed if it ends 
% with ".gz" or ".bz2"
%
%
% Alexander Barth


function fid=gzfopen(fname,permisson,machineformat)

global GZ_FILE_IDS GZ_FILE_NAMES

if isempty(GZ_FILE_NAMES)
  GZ_FILE_NAMES={};
end

[fname] = gread_tilde_expand(fname);


[tmp] = decompress(fname);
zipped = ~strcmp(tmp,fname);

if zipped
  fid = fopen(tmp,permisson,machineformat);

  GZ_FILE_IDS(end+1) = fid;
  GZ_FILE_NAMES{end+1} = tmp;
else
  fid = fopen(fname,permisson,machineformat);
end

