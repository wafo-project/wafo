function [y, tol1] = tay81cdf(x,a,b,tol)
%TAY81CDF Tayfun (1981) CDF of  breaking limited wave heights 
%                  /inf 
%   pdf = x/a^2 * | u* J_o(u/b)^(b^2)*J_o(xu) du *(1 - 4/pi*acos(b/x)*H(1-b/x) )
%                 /0
% where J_o and H(.) is the zero order Bessel function and the heavyside function
%
%  CALL:  [F tol]= tay81cdf(x,A,B,reltol) 
%
%            F = CDF
%            x = wave heights  0 <= x <= sqrt(2)*B
%            A = scale parameter
%            B = maximum limiting value (B->inf => pdf->Rayleigh )
%       reltol = relative tolerance (default 1e-3)
%           tol = absolute tolerance abs(int-intold)
%
% Tayfun set the scale parameter and  breaking limiting value to
%     A = Hrms = 2*sqrt(2*m0)*sqrt(1-0.734*eps2^4)
%     B = pi*tanh(k0*h)/(7*k0*2*sqrt(m0))
% where
%    k0 = the mean apparent wave number
%    mi = i'th spectral moment (i=0 => variance of the process)
%  eps2 = sqrt(m0*m2/m1^2-1) spectral bandwidth
%     h = water depth
%
%
%The size of F is the common size of X, A and B.  A scalar input   
%    functions as a constant matrix of the same size as the other input.
%
% Example:
% %  Hm0  = 7; Tp = 10; h = 150;
% %  S    = jonswap([],[Hm0,Tp]);
% %  eps2 = spec2char(S,'eps2');
% %  Tm02 = spec2char(S,'Tm02');
% %  a = Hm0/sqrt(2)*sqrt(1-0.734*eps2^4);
% %  k0   = w2k(2*pi/Tm02,[],h);
% %  b = pi*tanh(k0*h)/(7*k0*Hm0/2);
% %  x = linspace(0,Hm0)';
% %  plot(1:5,tay81cdf(1:5,a,b),x,cdfray(x,Hm0/2),'r')
%
% See also  tay81pdf, braylcdf

%   Reference:
%  M. Aziz Tayfun (1981) "Breaking limited waveheights"
% Journal of the Waterway, port and coastal and ocean division 
%      vol 107 No.2 pp. 59-69

%tested on:
%history
%
%  Per A. Brodtkorb 01.02.99

% TODO % not tested
  
if (nargin <4)|| isempty(tol), 
  tol=1e-2;%sqrt(eps);%relative accuracy of the estimates
end
if nargin <  3, 
    error('Requires at least three input arguments.'); 
end

[icode x a b] = iscomnsize(x,a,b.^2);

if ~icode 
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y=zeros(size(x));
tol1=y;
% Return NaN if A or B is not positive.
k1 = find(a <= 0 | b<= 0);
if any(k1) 
    tmp   = NaN;
    y(k1) = tmp(ones(size(k1)));
end

trace=false; %logical(0);%[];%used for debugging

k=find(a > 0 & x > 0 &x<=a.*sqrt(2*b) & b<inf&b> 0);
if any(k),
  %yk = y(k);
  xk = x(k); ak = a(k); bk = b(k);
  %yk = yk(:); 
  xk = xk(:); ak = ak(:); bk = bk(:);
  if 0,all((ak==ak(1)).*(bk==bk(1))),
%     [xk ind]=sort(xk);
%     gn=2;
%     [y(k), tol(k)]= gaussq('tay81pdf',[0;xk(1:end-1)] ,xk,tol,trace,gn,ak,bk,tol);
%     %[y(k), tol(k)]= asimpson('tay81pdf',0,xk,tol,trace,ak,bk);
%     y(k)=cumsum(y(k));
%     y(k)=y(k(ind)); %make sure the ordering is right
  else
    gn=2;
    [y(k), tol1(k)]= gaussq(@tay81pdf,0,xk,tol,[trace,gn],ak,bk,tol);
    %[y(k), tol(k)]= asimpson('pdfraymod',0,xk,tol,trace,ak,bk);
  end
end
if 0
  % This is a trick to get the html documentation correct.
  k = tay81pdf(1,1);
end

k2=find((a > 0 & y>1) | (x>=a.*sqrt(2*b) )& b<inf & b> 0);
if any(k2)
     y(k2)= ones(size(k2));
end

k4=find(a > 0 & x >= 0 & b==inf);
if any(k4), % special case b==inf -> rayleigh distribution
 y(k4) = cdfray(x(k4),a(k4)/sqrt(2));
end

% function z=traylfunzeros(x,a,b,n)
% % internal function which calculates the approximate zerocrossings of 
% %  the integrand:  x/a* J_o(u/sqrt(b))^b*J_1(x/a* u)
% nx=1:n;
% z1=a./x; 
% z0=sqrt(b);
% %k=find(mod(b,2)==0); % remove
% %z0(k)=pi/2*(n+1)*max(z1);
% z0=z0(:,ones(1,n)).* pi.*(nx(ones(size(x)),:)+0.25);
% z1=z1(:,ones(1,n)).* pi.*(nx(ones(size(x)),:)-0.25);
% z=sort([z0,z1],2);
% return


