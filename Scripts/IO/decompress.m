function [filename] = decompress(gzfile)


zipped = 0;

if length(gzfile) > 3
  if strcmp(gzfile(end-2:end),'.gz')
     zipped = 1;
  end
end

if length(gzfile) > 4
  if strcmp(gzfile(end-3:end),'.bz2')
     zipped = 2;
  end
end

filename = gzfile;

if zipped > 0
  tmpdir = getenv('TMPDIR');

  if isempty(tmpdir)
    tmpdir = '/tmp';
  end
  
  filename = tempname(tmpdir);

  if zipped == 1
    cmd = ['gunzip --stdout "' gzfile '"  >  "' filename '"'];
  else
    cmd = ['bunzip2 --stdout "' gzfile '"  >  "' filename '"'];
  end
  
  status = system(cmd);    
  if status ~= 0
    error(['decompress: command "' s '" failed ']);
  end
end
