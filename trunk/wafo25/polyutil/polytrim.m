function r = polytrim(p)
%POLYTRIM Trim polynomial by stripping off leading zeros.
%
%   R = POLYTRIM(P) trims the polynomial P by stripping off unnecessary
%   leading zeros.
%
%   P is a vector of coefficients in decreasing order.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:31:39
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   k = find(p);         % Find non-zero coefficients.
   if isempty(k)        % If non were found...
      r = 0;            % ...return zero polynomial...
   else                 % or else...
      k = min(k);       % ...get index of first...
      n = length(p);    % ...get length of vector...
      r = p(k:n);       % ...and assign output polynomial.
   end
