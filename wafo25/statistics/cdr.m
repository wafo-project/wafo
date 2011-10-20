function [I0,R] = cdr(I,varargin)
%CDR Complete Defining Relation
%
% CALL:  I0 = cdr(I);
%
%  I0 = Matrix of defining relations including all words.
%  R  = Integer defining the resolution.
%  I  = Matrix of generators.
%
% CDR uses the matrix of P generators to compute the complete defining
% relation for the design. This is useful to determine the confounding
% patterns of a two-level fractional design. The resolution of the design
% is also identified as the length of the shortest word in the defining
% relation. Any P words of I0 may be used as a genrator for the design.
% 
% Example:
%   I = sudg(8,4);
%   D = ffd(8,I);    % Fractional design for 8 variables in 2^(8-4) runs.
%   [I0,r] = cdr(I); % with a resolution IV.
%
% See also  sudg


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




error(nargchk(1,2,nargin))

if ischar(I),
 I = cl2cnr(I); % Transform from a column label into column nr.
end

I1 = abs(I);
p = size(I1,1);
k = max(I1(:));
if k-p<2
  warning('WAFO:CDR','This is not a design generator!')
  I0 = I;
elseif p>1,
  m = 2^p-1;
  n = size(I1,2);
  
  sgn = sign(I);
  sgn(sgn==0)=1;
  sgn = prod(sgn,2);
  sgn0 = ones(m,1);
  
  I0 = zeros(m,k);
  iz = 0;
  for ix = 1:p,
    iz      = iz+1;
    I0(iz,k+1-n:k) = sort(I1(ix,:));
    sgn0(iz) = sgn(ix);
    
    iz0     = iz;
    sgn0(iz0+1:2*iz0-1) = sgn0(iz0)*sgn0(1:iz0-1); 
    for iy = 1:iz0-1,
      iz = iz+1;
     % find values that are not in the intersection of I0(iz0,:) and I0(iy,:)
      tmp = setxor(I0(iz0,:),I0(iy,:));
      I0(iz,k+1-length(tmp):k)   = tmp;
    end
  end
  I0(:,k)= I0(:,k).*sgn0;
else
  I0 = I;
end
if 1,
  % Remove leading zeros
  [ix,iy] = find(I0~=0); 
  if ~isempty(iy)
    ix = min(iy);
    I0 = I0(:,ix:end);
    k = size(I0,2);
  end
else
  for ix=1:k,
    if any(I0(:,ix)~=0), I0 = I0(:,ix:k); k = k-ix+1; break,end
  end
end

if nargout>1,
  % Length of the shortest word in the defining relation 
  % is equal to the resolution.
  R = size(I0,2);
  [ix,iy] = find(I0==0);
  if any(iy),
    R = R-max(iy(:));   
  end
end

I0 = sortrows(I0);


k1 = find(I0==0);
if any(k1),
  % Sort so that zeros comes last on each row.
  I0(k1) = NaN;
  I0 = sort(I0,2);
  I0(isnan(I0))=0;
end

if nargin<2 && k<=50, % secret option
 I0 = cnr2cl(I0); % Transform into a column label
end
return



