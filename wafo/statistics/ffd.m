function [D,I] = ffd(n,I)
%FFD Two-level Fractional Factorial Design
%
% CALL:  [D,I] = ffd(n,I);
%        [D,I] = ffd(n,p);
%
%  D  = design matrix  size  2^(n-p) x n
%  n  = number of variables
%  I  = matrix of design generators. size p X q (default [])
%  p  = Number of generators
%
% FFD constructs a fractional factorial design for N variables in 
% 2^(N-P) runs. The matrix I contains P generators that define the
% design, see SUDG. If I is empty a full factorial design is returned.
% In general, increase in the degree of fractionation lowers the
% resolution of the best fraction and increases confounding between
% effects of various orders.  
%
% Examples
%   I1 = sudg(15,11); 
%   D1 = ffd(15,I1);  % Saturated Resolution III design.
% 
% See also  sudg

% Reference 
% Box, G.E.P, Hunter, W.G. and Hunter, J.S. (1978)
% Statistics for experimenters, John Wiley & Sons, pp 410

% Tested on: Matlab 5.3
% History:
% By Per A. Brodtkorb 16.03.2001

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


%error(nargchk(1,2,nargin))
narginchk(1,2)
if nargin<2||isempty(I); I = zeros(0,1); end
if length(I) ==1, I = sudg(n,I); end
if ischar(I),  I = cl2cnr(I); end


p    = size(I,1);
nmp  = n-p;
rows = 2^(nmp);
x1   = (0:rows-1)';
x2   = 1:nmp; %(nmp):-1:1;

% Construct a full factorial design for n-p variables
% in standard order.
D = bitget(x1(:,ones(1,nmp)),x2(ones(rows,1),:));

if p>0,
  II = abs(I);
  if any(II(:)>n), 
    error('Integers of the matrix of generators must be less or equal to %.0f ',n)
  end
  
  [ix,iy] = find(nmp<II & II <=n);
  if isempty(ix),
    ind = nmp:n;
  elseif length(ix)==p && all(diff(sort(ix(:))))
    ind     = I(sub2ind(size(I),ix,iy));
    [iz,iu] = sort(ix);
    ind     = ind(iu);
  else
    error('Something wrong with the matrix of generators, I')
  end
  D2 = D;
  D2((D2==0))=-1;
  D = [D zeros(rows,p)];
  for ix = 1:p,
    k = find(0<II(ix,:) & II(ix,:)<= nmp);
    if any(k),
      sgn = prod(sign(I(ix,k)))*ind(ix);
      D(:,abs(ind(ix))) = (prod(D2(:,II(ix,k)),2)*sgn>0);
    else
      error('Matrix of generators can not contain rows with only zeros')
    end
  end
end
if nargout>1,
   I = cnr2cl(I);
end
return

