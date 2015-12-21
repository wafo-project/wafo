%RNDDISCRETE Random sample
% 
% CALL  rnddiscrete(n, v, p)
%       rnddiscrete(v, p, r, c)
%       rnddiscrete(v, p, sz)
%
% Generate a row vector containing a random sample of size n from
% the univariate distribution which assumes the values in v with
% probabilities p. n must be a scalar.
%
% If r and c are given create a matrix with r rows and
% c columns. Or if sz is a vector, create a matrix of size
% sz.
% 
% See also rndempirical, pdfdiscrete, cdfdiscrete, invdiscrete

% Copyright (C) 1996, 1997  Kurt Hornik
%
% This file is part of WAFO.
%
% WAFO is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% WAFO is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with WAFO; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.


% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% pab 2007
% -translated from octave
% Description: Random deviates from a discrete distribution

function rnd = rnddiscrete(v, p, r, c)


tag1 = 'WAFO:RNDDISCRETE';
  if (nargin == 4)
    if (~(isscalar (r) && (r > 0) && (r == round (r))))
      error (tag1,'r must be a positive integer');
    end %if
    if ( ~(isscalar (c) && (c > 0) && (c == round (c))))
      error (tag1,'c must be a positive integer');
    end % if
    sz = [r, c];
  elseif (nargin == 3)
    % A potential problem happens here if all args are scalar, as
    % we can't distiguish between the command syntax. Thankfully this
    % case doesn't make much sense. So we assume the first syntax
    % if the first arg is scalar

    if (isscalar (v))
      sz = [1, floor(v)];
      v = p;
      p = r;
    else
      if (isscalar (r) && (r > 0))
        sz = [r, r];
      elseif (isvector(r) && all (r > 0))
        sz = r(:)';
      else
        error (tag1,'r must be a postive integer or vector');
      end %if
    end %if
  else
    error(tag1,'Too few inputs')
  end %if

  if (~ isvector (v))
    error (tag1,'v must be a vector');
  elseif (~isvector (p) || (length (p) ~= length (v)))
    error (tag1,'p must be a vector with length (v) elements');
  elseif (~ (all (p >= 0) && any (p)))
    error (tag1,'p must be a nonzero, nonnegative vector');
  end %if
  %I = ceil(n*rand(nn,1));
  np = numel(p);
  I = interp1(cumsum (p) / sum(p),1:np,rand(sz),'nearest','extrap');
  rnd = v (I); 
end %function
