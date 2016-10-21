%CDFEMPIRICAL Empirical cdf
%
%  CALL  F = cdfempirical(x,data)
% 
% For each element of x, compute the cumulative distribution
% function (CDF) at x of the empirical distribution obtained from 
% the univariate sample @var{data}.
%
% See also cdfdiscrete, pdfempirical, invempirical, rndempirical


% Copyright (C) 1996, 1997  Kurt Hornik
%
% CDFEMPIRICAL is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% CDFEMPIRICAL is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CDFEMPIRICAL; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% revised pab 2007
% -translated from octave to matlab
% -renamed from empirical_cdf to cdfemperical.
% Description: CDF of the empirical distribution

function cdf = cdfempirical (x, data)

  if (~isvector (data))
    error ('WAFO:CDFEMPIRICA','data must be a vector');
  end %if

  cdf = cdfdiscrete(x, data, ones (size (data)) / length (data));

end %function
