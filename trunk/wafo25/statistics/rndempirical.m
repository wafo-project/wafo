%RNDEMPIRICAL Bootstrap sample
% 
% CALL  rndempirical(n, data)
%       rndempirical(data, r, c)
%       rndempirical(data, sz)
%
% Generate a bootstrap sample of size @var{n} from the empirical
% distribution obtained from the univariate sample @var{data}.
%
% If r and c are given create a matrix with r rows and
% c columns. Or if sz is a vector, create a matrix of size
% sz.
%
% See also pdfempirical, cdfempirical, invempirical

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
% Description: Bootstrap samples from the empirical distribution

function rnd = rndempirical(data, r, c)

tag1 = 'WAFO:RNDEMPIRICAL';
  if (nargin == 2)
    if (isscalar(data))
      c = data;
      data = r;
      r = 1;
    end %if
  elseif (nargin ~= 3)
    error(tag1,'Too few input arguments')
  end %if

  if (~isvector (data))
    error (tag1,'data must be a vector');
  end %if

  rnd = rnddiscrete(data, ones (size (data)) / length (data), r, c);

end %function
