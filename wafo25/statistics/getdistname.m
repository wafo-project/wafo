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

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.



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
