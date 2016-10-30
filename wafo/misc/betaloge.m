function y = betaloge(z,w)
%BETALOGE Natural Logarithm of beta function.
%
% CALL y = betaloge(z,w)
%
%   
%     BETALOGE computes the natural logarithm of the beta
%     function for corresponding elements of Z and W.   The arrays Z and
%     W must be real and nonnegative. Both arrays must be the same size, 
%     or either can be scalar.  BETALOGE is defined as:
%  
%       y = LOG(BETA(Z,W)) = gammaln(Z)+gammaln(W)-gammaln(Z+W) 
%  
%     and is obtained without computing BETA(Z,W). Since the beta
%     function can range over very large or very small values, its
%     logarithm is sometimes more useful.
%     This implementation is more accurate than the BETALN implementation
%     for large arguments
%
% Example
%
%
%
% See also betaln, beta

% y = gammaln(z)+gammaln(w)-gammaln(z+w);
zpw = z+w;
y = stirlerr(z) + stirlerr(w) +0.5*log(2*pi)+(w-0.5).*log(w)+(z-0.5).*log(z)...
-stirlerr(zpw)-(zpw-0.5).*log(zpw);


% stirlings approximation:
%  (-(zpw-0.5).*log(zpw) +(w-0.5).*log(w)+(z-0.5).*log(z) +0.5*log(2*pi));

%!test assert(betaloge(1,1), 2.22044604925031e-016, 1e-12)
%!test assert(betaloge(20,1), -2.99573227355399, 1e-14)
%!test assert(betaloge(60,1), -4.09434456222215, 1e-14)
%!test assert(betaloge(100,1), -4.60517018598807, 1e-14)
%!test assert(betaloge(1000,1), -6.90775527898313, 1e-14)
