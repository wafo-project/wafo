function info = disectmfile(fname)
% DISECTMFILE Return functions called and scope of internal functions
%
% CALL info = disectmfile(file)
%
%  info = struct array with fieldnames:
%     .tag   : Indicating if a function is called or scope of internal function
%              The tag is a letter followed by an integer defining the indentation level.
%              The following letters means:
%           'M' = start scope of the main function
%           'N' = start scope of the nested function
%           'S' = start scope of the subfunction
%           'E' = end   scope of the internal function
%           'U' = external/internal function called
%     .line  : line   number of the TAG
%     .column: column number of the TAG
%     .name  : function name
%
% Example
%  if ismatlab(),
%    t = disectmfile('disectmfile')
%  end
%
% See also scan, mlint

% Tested on: matlab 7
% By pab 2007

str = mlint('-string','-calls',fname);
%str = evalc(sprintf('mlint(''-calls'',''%s'')', fname));
%names = regexp(str,'\d+\s*([AZa-z][A-Za-z0-9_]*)','tokens');
%info.names  = cellfun(@(x) x{1}, names, 'uniformoutput', false);
% if 0,
%   [info.tag,info.line,info.column,info.names] = strread(str,'%s %d %d %s');
% else
  pattern = '(?<tag>\w+)\s+(?<line>\d+)\s+(?<column>\d+)\s+(?<name>\w+)';
  info = regexp(str,pattern,'names');
% end