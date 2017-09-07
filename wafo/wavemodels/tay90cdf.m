function [y, tol1] = tay90cdf(x,a,b,tol)
%TAY90CDF Tayfun (1990) CDF of large wave heights 
%        /H             /1
%   F = |.5*H^3/A^2 *  |  z/sqrt(1-z)* exp(-(H/A)^2*(2-z))*Io((H/A)^2*z*b/2) dz dH
%       /0             /0
%                    
%   where Io(.) is the modified bessel funcition of zero order 
%
%  CALL:  [F, tol]= tay90cdf(H,A,B,reltol);
%
%            F = CDF
%            H = wave height
%            A = scale parameter (=sqrt(m0))
%            B = correlation between At and Ac (B->1 => pdf->Rayleigh )
%       reltol = relative tolerance (default 1e-3)
%          tol = absolute tolerance abs(int-intold)
%
% The size of F is the common size of X, A and B.  A scalar input   
% functions as a constant matrix of the same size as the other input.
%
%  Tayfun calculates the correlation from the Spectral density as:
%
%  B = sqrt( [int S(w)*cos(w*T) dw]^2+[int S(w)*sin(w*T) dw]^2 )/m0
%  
%  where T = Tm02, mean zero crossing period.
%  Note that this is the same as the groupiness factor, Ka, found in
%  spec2char
%
% Example:
%  S = jonswap;
%  a = sqrt(spec2mom(S,1));
%  b = spec2char(S,'Ka');
%  h = linspace(0,7*a)';
%  plot(h,tay90cdf(h,a,b))
%
% See also  tay90pdf, spec2char, gaussq


%   References:
% [1]  M. Aziz Tayfun (1990) "Distribution of large waveheights"
% Journal of the Waterway, port and coastal and ocean division 
%      vol 116 No.6 pp. 686-707
% [2]  M. Aziz Tayfun (1981) "Distribution of crest-to-trough waveheights"
% Journal of the Waterway, port and coastal and ocean division 
%      vol 107 No.3 pp. 149-158

%tested
% history:
% revised pab 01.12.2000
%  -fixed a bug in the normalization
%  Per A. Brodtkorb 21.02.99

%error(nargchk(3,4,nargin))
narginchk(3,4)
if (nargin <4)|| isempty(tol), 
  tol=1e-3; %relative accuracy of the estimates % sqrt(eps);%
end

[icode, x, a, b] =iscomnsize(x,a,b);

if ~icode 
     error('h, a and b must be of common size or scalar.');
end
xsiz = size(x); % remember old size
x=x(:); a=a(:);b=b(:);

% Initialize Y to zero.
y=zeros(size(x));
tol1=y;
% Return NaN if A is not positive.
k1 = find(a <= 0| abs(b)>1);
if any(k1) 
    tmp   = NaN;
    y(k1) = tmp(ones(size(k1)));
end



gn    = 2;% # of of base points to start the integration with
trace = 0;%used for debugging


if 0
  % This is a trick to get the html documentation correct.
  k = tay90pdf(1,1,0.1);
end

k=find(a > 0 & x >0 & abs(b)<1);
if any(k),
  %yk=y(k);
  xk = x(k); ak = a(k); bk = b(k);
  
 [yk, tol1(k)]= gaussq(@tay90pdf,0,xk,  tol,[trace gn],ak,bk,tol);
 % yk=yk./ak;
  
  k2=find(yk>1 );
   if any(k2)
     yk(k2)= ones(size(k2));
  end
 y(k)=yk;
end

k4=find(a > 0 & x >= 0 & abs(b)==1);
if any(k4), % special case b==1 -> rayleigh distribution
 y(k4) = cdfray(x(k4),2*a(k4));
end 
y=reshape(y,xsiz);

