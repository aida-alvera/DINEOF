% gzfclose(fid)
% 
% same as fclose except that a decompressed temporary file 
% is deleted if the file ends with ".gz" or ".bz2"
%
%
% Alexander Barth, 2008-03-19

function fid = gzfclose(fid)

global GZ_FILE_IDS GZ_FILE_NAMES

fid = fclose(fid);

index = [];
if ~isempty(GZ_FILE_IDS)
  index = find(GZ_FILE_IDS == fid);
end

if ~isempty(index)
  tmp = GZ_FILE_NAMES{index};
  delete(tmp);

  GZ_FILE_IDS(index)=[];
%  GZ_FILE_NAMES={GZ_FILE_NAMES{1:index-1}  GZ_FILE_NAMES{index+1:end}};
  
% for compatability witch octave  
  GZ_FILE_NAMES(index:end-1) = GZ_FILE_NAMES(index+1:end);
  GZ_FILE_NAMES = GZ_FILE_NAMES(1:end-1)  ;
end

