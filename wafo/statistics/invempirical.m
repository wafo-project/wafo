% INVEMPIRICAL Empirical quantile. 
% 
% CALL x = invempirical(F, data)
%
% For each element of F, compute the quantile (the inverse of the
% CDF) at F of the empirical distribution obtained from the
% univariate sample data.
%
% See also invdiscrete, pdfempirical, cdfempirical, rndempirical


% Copyright (C) 1996, 1997  Kurt Hornik
%
% INVEMPIRICAL is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% INVEMPIRICAL is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with INVEMPIRICAL; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% pab 2007
% -translated from octave
% Description: Quantile function of the empirical distribution

function inv = invempirical (x, data)

  if (~isvector (data))
    error ('WAFO:INVEMPIRICAL','data must be a vector');
  end %if

  inv = invdiscrete(x, data, ones (size (data)) / length (data));

  end %function
