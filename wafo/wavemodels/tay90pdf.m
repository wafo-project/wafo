function [y, tol1] = tay90pdf(x,a,b,tol)
%TAY90PDF Tayfun (1990) PDF of large wave heights 
%                      /1
%   f = .5*H^3/A^2 *  |  z/sqrt(1-z)* exp(-(H/A)^2*(2-z))*Io(  (H/A)^2*z*b/2) dz
%                     /0
%                    
%   where Io(.) is the modified bessel funcition of zero order 
%
%  CALL:  [f tol]= tay90pdf(H,A,B,reltol)
%        
%            f = pdf
%            H = wave height
%            A = scale parameter (=sqrt(m0))
%            B = correlation between trough and crest amplitude squared
%                i.e. Cov(At^2, Ac^2) approx -R(T/2)  (B->1 => pdf->Rayleigh )      
%       reltol = relative tolerance default=1e-3
%          tol = absolute tolerance abs(int-intold)
%
% The size of f is the common size of X, A and B.  A scalar input   
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
%  plot(h,tay90pdf(h,a,b))
%
% See also  tay90cdf, spec2char, gaussq



%   Reference:
% [1]  M. Aziz Tayfun (1990) "Distribution of large waveheights"
% Journal of the Waterway, port and coastal and ocean division 
%      vol 116 No.6 pp. 686-707
% [2]  M. Aziz Tayfun (1981) "Distribution of crest-to-trough waveheights"
% Journal of the Waterway, port and coastal and ocean division 
%      vol 107 No.3 pp. 149-158

%tested on: 
%history:
% revised pab 26.07.2001
% - added forgotten ) in error(nargchk(3,4,nargin))
% revised IR, JR, PAB
%  raylpdf -> wraylpdf
% revised pab 07.03.2000
% - updated error in help header
% revised pab 28.02.2000
%  - fixed some bugs
%  Per A. Brodtkorb 21.02.99


%error(nargchk(3,4,nargin))
narginchk(3,4)
if (nargin <4)|| isempty(tol), 
  tol=1e-3; %relative accuracy of the estimates % sqrt(eps);%
end

[icode x a b] = iscomnsize(x,a,b);

if ~icode 
  error('h, a and b must be of common size or scalar.');
end
xsz = size(x);
x=x(:); a=a(:);b=b(:);
%a
%b
% Initialize Y to zero.
y=zeros(size(x));
tol1=y;
% Return NaN if A is not positive.
k1 = find(a <= 0| abs(b)>1);
if any(k1) 
    tmp   = NaN;
    y(k1) = tmp(ones(size(k1)));
end

if 0
  % This is a trick to get the html documentation correct.
  k = tay90fun(1,1);
end

gn=2;% # of of base points to start the integration with
trace=0;%used for debugging
%size(x)
%size(a)
%size(y)
k=find(a > 0 & x >0 & abs(b)<1);
%size(k)
if any(k),
  yk=y(k);xk = x(k); ak = a(k);
  
 [yk, tol1(k)]= gaussq(@tay90fun,0,ones(size(yk)),  [tol 4],[trace gn],xk./ak,b(k));
  yk=yk./ak;
  
  k2=find(yk<0 );
   if any(k2)
     disp('Some less than zero')
     yk(k2)= zeros(size(k2));
  end
 y(k)=yk;
end

k4=find(a > 0 & x >= 0 & abs(b)==1);
if any(k4), % special case b==1 -> rayleigh distribution
 y(k4) = pdfray(x(k4),2*a(k4));
end 
y=reshape(y,xsz);






