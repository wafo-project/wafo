
% SCAN Scans and extracts function names from an m file
% 
%  CALL names = scan(file,patterns)
%
%
% Arguments
%   file        -   full path to an m file to scan
%   patterns    -   A cell array of regular expressions for function names
%                   to extract
%
% SCAN Scans and extracts function names from an m file
% using regular expression patterns.
%
%
function names = scan(fname, patterns)
  if nargin<2
    patterns = '';
  end
  if iscell(patterns)
    patterns = sprintf('(%s)|',patterns{:});
  end
  str = mlint('-string','-calls',fname);
  %str = evalc(sprintf('mlint(''-calls'',''%s'')', fname));
  names = regexp(str,'\d+\s*([AZa-z][A-Za-z0-9_]*)','tokens');
  names = cellfun(@(x) x{1}, names, 'uniformoutput', false);
  if ~isempty(patterns)
    i = cellfun(@(x)~isempty(regexp(x, patterns, 'once' )), names);
    names = names(i);
  end
end
