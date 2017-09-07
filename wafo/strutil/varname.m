function str = varname(varargin)
%VARNAME Return the name of the variable as a string
%
% CALL: str = varname(x)
%
% Example:
% x2 = 'test';
% assert(varname(x2), {'x2'})
% assert(varname(1), {''})
  
  %error(nargchk(1,inf,nargin));
  narginchk(1,inf)
  str = cell(1,nargin);
  for ix=1:nargin,
    str{ix} = inputname(ix);
  end
  return