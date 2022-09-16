% files = filelist(pattern);
% create a list of a file names matching the pattern (include wild-cards such as *)
% expands ~ to home directory and sorts all files alphabetically


function files = filelist(pattern);

pattern = gread_tilde_expand(pattern);

i = regexp(pattern,'[#(]');

if ~isempty(i)
  subset = pattern(i:end);
  pattern = pattern(1:i-1);
else
  subset = '';
end

%[basedir] = fileparts (pattern);
%d = dir(pattern);
l = uglob(pattern);

files = cell(1,length(l));
for i=1:length(l)
  files{i} = [l{i} subset];
end

files = sort(files);


