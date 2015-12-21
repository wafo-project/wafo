% PDFEMPIRICAL Empirical pdf
%
% CALL f = pdfempirical(x, data)
%
% For each element of x, compute the probability density function
% (PDF) at x of the empirical distribution obtained from the
% univariate sample DATA.
%
% See also pdfdiscrete, cdfempirical, invempirical, rndempirical

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
% Description: PDF of the empirical distribution

function pdf = pdfempirical(x, data)

tag1 = 'WAFO:PDFEMPIRICAL';
  if (~ isvector (data))
    error (tag1,'data must be a vector');
  end %if

  pdf = pdfdiscrete(x, data, ones (size (data)) / length (data));

end %function
