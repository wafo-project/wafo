function R = var2corr(S)
%VAR2CORR Variance matrix to correlation matrix conversion.
%
%   R = VAR2CORR(S) returns the correlation matrix, which is the variance
%   matrix normalized so that the diagonal elements are ones or zeros and
%   the off-diagonal elements are correlations:
%
%      R(i,j) = S(i,j) / sqrt( S(i,i) * S(j,j) )
%
%   R will always have ones on the diagonal.  S must be positively
%   semidefinite.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-06-19 12:40:32 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));        % check number of input args

   [m, n] = size(S);
   if m ~= n
      error('Input matrix must be square.');
   end

   R = zeros(size(S));                  % initialize output matrix
   v = diag(S);                         % variances
   k = find(v);                         % find non-zeros on diagonal

   if ~isempty(k)
      S_sub = S(k,k);                   % get non-singular submatrix
      D = diag(sqrt(v(k)));             % standard deviations
      R(k,k) = D\S_sub/D;               % correlation submatrix
   end
