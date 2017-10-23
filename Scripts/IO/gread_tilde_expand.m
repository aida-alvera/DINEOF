% replace ~ by home dir 
% (function build-in in octave)
%
% Alexander Barth, 2008-07-16

function [fname] = gread_tilde_expand(fname)

homedir = getenv('HOME');
if ~isempty(homedir) & fname(1) == '~'
  fname = strrep(fname,'~',homedir);
end
