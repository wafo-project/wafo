function [varargout]=deriv(t,k)
%DERIV  4th, 6th, 8th and 10th derivatives of the kernel function.
%
% CALL:  [y4,y6,y8,y10] = deriv(t,kernel)
%
%   kernel = 'gaussian'      - Gaussian kernel (Currently the only
%                              supported kernel)
%
%   DERIV(T,K) finds the derivatives of the kernel K 
%              at the point T. 
%
% Example: [y4,y6,y8,y10]=deriv(0:3);

%tested on: matlab 5.3
%History
%revised pab Aug2005
% made it general for derivatives of any order
%revised pab dec2003
% added todo  
%revised pab 16.10.1999
%  updated to matlab 5.x + documentation
% by             Christian C. Beardah 1995

  
% TODO % Add support for other kernels than the Gaussian
  
if nargin<2||isempty(k)
  k='gauss';
end

switch lower(k(1:4))
  case {'gaus'},
    phi0 = exp(-0.5*t.^2)/sqrt(2*pi);
    %if 1
      % New call by pab Aug 2005
      p4 = [1 0 -6 0 +3];
      varargout{1} = polyval(p4,t).*phi0;

      pn = p4;
      for ix = 2:nargout
        pnp1 = polyadd([-pn,0],polyder(pn));
        pnp2 = polyadd([-pnp1,0],polyder(pnp1));
        varargout{ix} = polyval(pnp2,t).*phi0;
        pn = pnp2;
      end
    
%     else % old call kept just in case 
%       y4=(t.^4-6*t.^2+3).*phi0;
%       if nargout>1,
%         y6=(t.^6-15*t.^4+45*t.^2-15).*phi0;
%       end
%       if nargout>2;
%         y8=(t.^8-28*t.^6+210*t.^4-420*t.^2+105).*phi0;
%       end;
%       if nargout>3,
%         y10=(t.^10-45*t.^8+630*t.^6-3150*t.^4+4725*t.^2-945).*phi0;
%       end;
%     end
otherwise,
  error('Kernel not suported')
end



function r = polyadd(p, q)
%POLYADD Add polynomials.
%
%   R = POLYADD(P, Q) adds the polynomials whose coefficients are the
%   elements of the vectors P and Q.
%
%   See also POLYADD, POLYSUB, POLYMUL, POLYDIV, POLYPOW.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:36:21
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(2, 2, nargin));

   m = length(p);
   n = length(q);
   l = max(m, n);

   r = zeros(1, l);                     % Initialize output.
   r(l-m+1:l) = p;                      % Insert first polynomial.
   r(l-n+1:l) = r(l-n+1:l) + q;         % Add second polynomial

   r = polytrim(r);                     % Trim the result.
   
   
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
