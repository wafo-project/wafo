function [ef, id]=yates(y,varargin)
%YATES Calculates main and interaction effects using Yates' algorithm.
%
% CALL:  [ef, id] = yates(y);
%
%  id = identification vector of main and interaction effects.
%  ef = vector of average response, main effects and interaction effects.
%  y  = calculated response from a two-level complete factorial design.
%
% YATES applies the  Yates' algorithm to the responses Y to obtain 
% average response, main effects and interaction effects. Y is assumed to
% be arranged in what is called standard order. (The order of the actual
% running should, of course, be random).  EF(1,:) is the
% average response and EF(2:end,:) contain the main effects and
% interaction effects corresponding to the vector ID.
% 
% YATES may also be used in analyzing data from any 2^(K-P) fractional
% factorial design. The algorithm is applied in the usual way to any ambedded
% complete factorial in K-P factors, i.e., the responses must be
% rearranged so that K-P factors is a complete factorial in standard order.
% Then associate the calculated effects with their appropriate aliases
% using ALIAS. 
%
% Example:
%   D = ffd(3);                    % complete 2^3 design in standard order.
%   y = [60 72 54 68 52 83 45 80]; % Responses to design D.
%   [ef,id] = yates(y);
%
%   I = sudg(7,4);
%   D1 = ffd(7,I);                  % 2^(7-4) design in standard order.
%   y1 = [69 52 60 83 71 50 59 88]; % Responses to design D1.
%   [ef1,id1] = yates(y1);
%   alias(cdr(I),3)                 % associate with id1 
%
% See also  fitmodel, getmodel, alias, cdr, sudg, ffd ,plotresponse, nplot


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



% Reference 
% Box, G.E.P, Hunter, W.G. and Hunter, J.S. (1978)
% Statistics for experimenters, John Wiley & Sons, pp 342

% Tested on: Matlab 5.3
% History:
% By Per A. Brodtkorb 16.03.2001

error(nargchk(1,2,nargin))
sz = size(y);
n  = length(y); 
if prod(sz) == n, 
  y = y(:);       % Make sure it is a column vector
else   
  n = sz(1);      % Number of runs
end

k = log2(n);      % Number of variables.
if round(k)~=k, error('The length of y must be in power of two'), end

% Yates algorithm:
%--------------------------
ef   = y;
ind2 = 2:2:n;
ind1 = 1:2:n-1;
for ix=1:k
  ef = [ef(ind2,:)+ef(ind1,:); ef(ind2,:)-ef(ind1,:)];
end
ef = ef*(2/n);
ef(1,:) = ef(1,:)/2;


if nargout>1,
  id = zeros(n-1,k);
  iz = 0;
  for ix = 1:k,
    iz       = iz+1;
    id(iz,1) = ix;
    iz0      = iz;
    for iy = 1:iz0-1,
      iz = iz+1;
      id(iz,:)   = id(iy,:);
      ind        = find(id(iy,:)==0, 1 );
      id(iz,ind) = ix;
    end
  end
  if nargin>1, % secret option 
    % Sort effects
    [id, ind] = sortrows(fliplr(id));
    id = fliplr(id);
    ef(2:end,:) = ef(ind+1,:);   
  end
  % String representation
  id = cnr2cl(id);
  %str0=[' ',char(65:90) char(97:122)]; % characters A - Z a-z
  %id = str0(id+1);
end


return




