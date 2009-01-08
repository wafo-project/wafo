function y = ochi98cdf(x,a,b)
%OCHI98CDF Ochi's (1998) CDF of peaks and troughs of non-gaussian processes 
%
%  CALL:  F = ochi98cdf(x,A1,A2) 
%
%           F = CDF
%           x = quantiles
%      A1, A2 = scale parameters for troughs and crests, respectively.
%  where
%      A1 = sqrt(2*s1^2/(1-(.4*(c1/s1)^2)
%      A2 = sqrt(2*s2^2/(1-(.4*(c2/s2)^2)
%      si = sigmai/(1+L2 
%
%  The size of F is the common size of X, A1 and A2.  A scalar input   
%  functions as a constant matrix of the same size as the other input.
%
% Example:
%  x =linspace(0,10).'; a1 = 2.5; a2 = 2;
%  plot(x,ochi98cdf(x,a1,a2),x,cdfray(x,sqrt(a1^2+a2^2))) 
% 
% See also  cdfnorm, cdfray



%  Reference:
%       [1]  Michel K. Ochi, (1998)
%       "Probability distributions of peaks and troughs of non-gaussian processes"
%        Probabilistic Engineering Mechanics, Vol 13, No 4, pp  291-298

% tested on:
% history:
% revised pab 29.02.2000
% changed name to ochi98cdf
%  Per A. Brodtkorb 17.02.99
  
% TODO % needs explanation to the scaling
if nargin <  3, 
    error('Requires at least three input arguments.'); 
end

[icode x a b] = iscomnsize(x,a,b);

if ~icode 
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y=zeros(size(x));

% Return NaN if A or B is not positive.
k1 = find(a <= 0|b<= 0);
if any(k1) 
    tmp   = NaN;
    y(k1) = tmp(ones(size(k1)));
end

k=find(a > 0 & b>0 & x > 0 & x < inf);
if any(k),
   r1 = a(k).^2; r2 = b(k).^2;xk = x(k);
 y(k)=1-1./(r1+r2).*(r1.*exp(-xk.^2./r1)+   r2.*exp(-xk.^2./r2)+...
           2*xk.*sqrt(pi*r1.*r2./(r1+r2)).*exp(-xk.^2./(r1+r2))...
    .*(cdfnorm(sqrt(2*r2./r1./(r1+r2)).*xk )    - cdfnorm(-sqrt(2*r1./r2./(r1+r2)).*xk ) ));
 end

k=find(a > 0 & (b==0 | a==b));
if any(k),
  y(k)=cdfray(x(k),a(k)/sqrt(2));
end
k=find(a == 0 & b>0);
if any(k),
  y(k)=cdfray(x(k),b(k)/sqrt(2));
end
k=find(b == 0 & a>0);
if any(k),
  y(k)=cdfray(x(k),a(k)/sqrt(2));
end
k=find(a > 0 & b>0 & x ==inf);
if any(k),
  y(k)=ones(size(y(k)));
end
