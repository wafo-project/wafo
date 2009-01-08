function name = getdistname(dist)
% GETDISTNAME Return the distribution name
%
% CALL name = getdistname(dist)
% 
%   dist = function handle or string to one of the cdfXXX,pdfXXX, prbXXX, 
%          invXXX and fitXXX functions
%
% Example
%  model = getdistname('pdfgev')
% 
% See also


switch class(dist)
  case {'char'}
    %dist = dist;
  case {'function_handle'} % OK
    dist = func2str(dist);
  otherwise
    error('Distribution is unspecified')
end

if numel(dist)>3 && (strncmpi(dist,'cdf',3) || strncmpi(dist,'pdf',3) || strncmpi(dist,'inv',3) ||...
    strncmpi(dist,'fit',3) || strncmpi(dist,'rnd',3) || strncmpi(dist,'prb',3))
  name = dist(4:end);
else
  name = dist; 
end
