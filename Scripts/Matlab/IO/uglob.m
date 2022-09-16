function f = uglob(pattern);

[status,output] = system(['ls -1 ' pattern]);
files = strsplit(output,sprintf('\n'));

f = {};

for i = 1:length(files)
  if length(files{i}) > 0
    f{end+1} = files{i};
  end
end
