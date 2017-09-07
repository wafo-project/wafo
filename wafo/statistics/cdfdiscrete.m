% CDFDISCRETE Discrete cdf
%
%  CALL F = cdfdiscrete(x,data, p)
% 
% For each element of x, compute the cumulative distribution
% function (CDF) at x of a univariate discrete distribution which
% assumes the values in data with probabilities P.
%
% Example
%  n = 50;        % Number of objects to test in a day
%  p = 0.2;      % Defect probability
%  defects = 0:n; % All possible outcomes of a day
%  f = pdfbin(defects,n,p);
%  F = cdfbin(defects,n,p);
%  xi = linspace(0,n,200);
%  fi = cdfdiscrete(xi,defects,f);
%  stairs(defects,F), hold on, 
%  plot(xi,fi,'r'),hold off, shg
% 
% See also cdfempirical, pdfdiscrete, invdiscrete, rnddiscrete



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

% -*- texinfo -*-
% @deftypefn {Function File} {} discrete_cdf (@var{x}, @var{v}, @var{p})
% For each element of @var{x}, compute the cumulative distribution
% function (CDF) at @var{x} of a univariate discrete distribution which
% assumes the values in @var{v} with probabilities @var{p}.
% @end deftypefn

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% pab 2007
% -translated from octave
% Description: CDF of a discrete distribution

function cdf = cdfdiscrete(x, v, p)

%error(nargchk(3,3,nargin))
narginchk(3,3)
%   if (nargin ~= 3)
%     print_usage ();
%   end %if

  sz = size (x);

  tag1 = 'WAFO:CDFDISCRETE';
  if (~isvector(v))
    error(tag1,'v must be a vector');
  elseif (~isvector(p) || (length (p) ~= length (v)))
    error (tag1,'p must be a vector with length (v) elements');
  elseif (~(all(p >= 0) && any (p)))
    error (tag1,'p must be a nonzero, nonnegative vector');
  end %if

  n = numel (x);
  m = length (v);
  x = reshape (x, n, 1);
  v = reshape (v, 1, m);
  p = reshape (p / sum (p), m, 1);

  cdf = zeros (sz);
  k = find (isnan (x));
  if (any (k))
    cdf (k) = NaN;
  end %if
  k = find (~isnan(x));
  if (any (k))
    n = length (k);
    cdf (k) = ((x(k) * ones (1, m)) >= (ones (n, 1) * v)) * p;
  end %if

  end %function
