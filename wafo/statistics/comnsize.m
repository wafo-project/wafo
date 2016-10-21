function [csize,varargout] = comnsize(varargin)
% COMNSIZE Calculates common size of all non-scalar arguments.
%
% CALL:  [csize,y1,y2,...,yN] = comnsize(x1,x2,...,xN);
%
%  csize = vector giving the common size of all non-scalar arguments.
%          If the non-scalar xi do not match in size then csize = nan.
%  y1,...,yN = Same as x1,...,xN, except that scalars are transformed to
%              a constant matrix with same size as the other inputs.
%  x1,...,xN = Input arguments.
%
%  COMNSIZE returns the common size of all non-scalar arguments and makes
%  sure that the output arguments Y1,Y2,...,YN, are of common size.
%  This is useful for implementing functions where arguments can either
%  be scalars or of common size. 
%
%  NOTE:  If the csize is nan, then yi = xi.
%
% Examples: 
%   A = rand(4,5);B = 2;C = rand(4,5);
%   [csiz,A1,B1,C1] = comnsize(A,B,C);
%   csiz = comnsize(A,1:2);

% Copyright (C) 1999, 2007 Per A. Brodtkorb
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

% Tested on: matlab 5.3,7.4
% History:
% revised pab aug 2007
% - vectorized further.
% revised pab 23.10.2000
%  - New name comnsize
%  - the input arguments can have a any dimension.
%  - Updated documentation
%  - Vectorized the calculations a bit more.
% Based on common_size.m by
%  KH <Kurt.Hornik@ci.tuwien.ac.at> Created 15 October 1994

Np   = nargin;
Nout = max(nargout-1,0);
if Nout>Np, 
  error('The number of output arguments is too large.')
end

vsz  = cellfun(@size, varargin,'UniformOutput', false);
dims = cellfun(@length,vsz);


is2D = (dims==2);
Ndim = max(dims);

sz           = ones(Np,Ndim);
sz(is2D,1:2) = cat(1,vsz{is2D});

for ix = find(~is2D)
  sz(ix,1:dims(ix)) = vsz{ix};
end

csize     = max(sz,[],1);        % find common size

% Check that non-scalar arguments match in size
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

isscalar   = all((sz == 1).'); % mask to scalars
iscomnsize = all((sz == csize(ones(Np, 1),:)).');
if any(~isscalar & ~iscomnsize)
  csize = nan;
end

if Nout>0
  nonscalar = ~isscalar(1:Nout); % Mask to nonscalar arguments in output

  if (any(isnan(csize)) || all(isscalar) || all(nonscalar))
    varargout(1:Nout) = varargin(1:Nout);
  else
    if any(nonscalar)
      varargout(nonscalar) = varargin(nonscalar);
    end
    % Make sure the scalar arguments are of size csiz
    for ix = find(isscalar),
      varargout{ix} = varargin{ix}(ones(csize));
    end
  end
end