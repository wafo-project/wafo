% INVDISCRETE Disrete quantile
%
% CALL x = invdiscrete(F, data,p)
%
% For each component of F, compute the quantile (the inverse of
% the CDF) at F of the univariate distribution which assumes the
% values in data with probabilities p.
%
% See also invempirical, cdfdiscrete, invdiscrete, rnddiscrete

% Copyright (C) 1996, 1997  Kurt Hornik
%
% INVDISCRETE is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% INVDISCRETE is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with INVDISCRETE; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% pab 2007
% -translated from octave
% Description: Quantile function of a discrete distribution

function inv = invdiscrete(x, v, p)

%error(nargchk(3,3,nargin))
narginchk(3,3)  
  sz = size (x);
  tag1= 'WAFO:INVDISCRETE';
  if (~isvector (v))
    error (tag1,'v must be a vector');
  elseif (~ isvector (p) || (length (p) ~= length (v)))
    error (tag1,'p must be a vector with length (v) elements');
  elseif (~ (all (p >= 0) && any (p)))
    error (tag1,'p must be a nonzero, nonnegative vector');
  end %if

  n = numel (x);
  x = reshape (x, 1, n);
  m = length (v);
  v = sort (v);
  s = reshape (cumsum (p / sum (p)), m, 1);

  % Allow storage allocated for P to be reclaimed.
  %p = [];

  inv = NaN * ones (sz);
  k = find (x == 0);
  if (any (k))
    inv(k) = -Inf;
  end %if
  k = find (x == 1);
  if (any (k))
    inv(k) = v(m) * ones (size (k));
  end %if
  k = find ((x > 0) & (x < 1));
  if (any (k))
    n = length (k);

    % The following loop is a space/time tradeoff in favor of space,
    % since the dataset may be large.
    %
    % Vectorized code is:
    %
    %     inv(k) = v(sum ((ones (m, 1) * x(k)) > (s * ones (1, n))) + 1);

    for q = 1:n
      inv(k(q)) = v(sum (x(k(q)) > s) + 1);
    end %for
  end %if

end %function


