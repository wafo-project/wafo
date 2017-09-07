% PDFDISCRETE Discrete PDF
% 
% CALL f =  pdfdiscrete(x, v, p)
%
% For each element of x, compute the probability density function
% (PDF) at x of a univariate discrete distribution which assumes
% the values in v with probabilities p.
%
% See also pdfempirical, cdfdiscrete, invdiscrete, rnddiscrete


% Copyright (C) 1996, 1997  Kurt Hornik
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% pab 2007
% -translated from octave
% Description: pDF of a discrete distribution

function pdf = pdfdiscrete(x, v, p)

%error(nargchk(3,3,nargi))
narginchk(3,3)
  sz = size (x);
tag1 = 'WAFO:PDFDISCRETE';
  if (~isvector (v))
    error (tag1,'v must be a vector');
  elseif (~isvector (p) || (length (p) ~= length (v)))
    error (tag1,'p must be a vector with length (v) elements');
  elseif (~(all (p >= 0) && any (p)))
    error (tag1,'p must be a nonzero, nonnegative vector');
  end %if

  n = numel (x);
  m = length (v);
  x = reshape (x, n, 1);
  v = reshape (v, 1, m);
  p = reshape (p / sum (p), m, 1);

  pdf = zeros (sz);
  k = find (isnan (x));
  if (any (k))
    pdf (k) = NaN;
  end %if
  k = find (~isnan (x));
  if (any (k))
    n = length (k);
    pdf (k) = ((x(k) * ones (1, m)) == (ones (n, 1) * v)) * p;
  end %if

end %function
