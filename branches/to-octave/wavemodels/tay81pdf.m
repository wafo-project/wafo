function [y, tol1] = tay81pdf(x,a,b,tol)
%TAY81PDF Tayfun (1981) PDF of  breaking limited wave heights 
%              /inf
%  f = x/a^2 * | u* J_o(u/b)^(b^2)*J_o(x/a* u) du *(1 - 4/pi*acos(b*a/x)*H(x/a-b) )
%              /0
% where J_o and H(.) is the zero order Bessel function and the heavyside function
%
%  CALL:  [f tol]= tay81pdf(H,A,B,reltol) 
%            f = pdf  
%            x = quantiles  0 <= x <= sqrt(2)*B
%            A = scale parameter
%            B = maximum limiting value (B->inf => pdf->Rayleigh )
%       reltol = relative tolerance (default 1e-3)
%          tol = absolute tolerance abs(int-intold)
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
%  The size of F is the common size of X, A and B.  A scalar input   
%  functions as a constant matrix of the same size as the other input.
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
% %  plot(x,tay81pdf(x,a,b),x,pdfray(x,Hm0/2),'r')
%
% See also  tay81cdf, braylpdf

%   Reference:
%  M. Aziz Tayfun (1981) "Breaking limited waveheights"
% Journal of the Waterway, port and coastal and ocean division 
%      vol 107 No.2 pp. 59-69

%tested on: 
%history:
%  Per A. Brodtkorb 01.02.99

% TODO % not tested
if (nargin <4)|| isempty(tol), 
  tol=1e-3; %relative accuracy of the estimates % sqrt(eps);%
end
if nargin <  3, 
    error('Requires at least three input arguments.'); 
end

[icode x a b] = iscomnsize(x,a,b.^2);

if ~icode 
    error('Requires non-scalar arguments to match in size.');
end
[N M]=size(x);
x=x(:); a=a(:);b=b(:);

% Initialize Y to zero.
y=zeros(size(x));
tol1=y;
% Return NaN if A or B is not positive.
k1 = find(a <= 0|b<= 0);
if any(k1) 
    tmp   = NaN;
    y(k1) = tmp(ones(size(k1)));
end


count_limit=20;
gn=2^3;% # of of base points to start the integration with
trace=false; %logical(0);%used for debugging

k=find(a > 0 & x >0 &x<a.*sqrt(2*b) & b<inf&b> 0);
if any(k),
  yk=y(k);xk = x(k); ak = a(k); bk = b(k);
  tmpy=yk;    tmpt=tmpy; %size(xk)
  uk_inf = traylfunzeros(xk,ak,bk,count_limit*10);%
  uk_inf = [zeros(size(k)) uk_inf(:,3:4:end)];
              
  if 0  
    uk_inf= uk_inf(:,3:4:end);
    k1=find( uk_inf(:,1)<=7);
    if any(k1),
      tmpy=uk_inf(:,1);
      uk_inf(k1,:)=uk_inf(k1,:)+7-tmpy(k1,ones(1,size(uk_inf,2)));
    end
    %uk_inf=.5*ak.*sqrt(2.*bk);
    %uk_inf(1:10,1:4), size(uk_inf),size(k)
    [yk, tol1(k)]= gaussq('tay81fun',0,uk_inf(k),  tol/2,trace,gn,xk./ak,bk);
  end
  if 0
    % This is a trick to get the html documentation correct.
    k = tay81fun(1,1,2);
  end
  
  
  % Break out of the iteration loop for three reasons:
  %  1) the last update is very small (compared to int)
  %  2) the last update is very small (compared to tol)
  %  3) There are more than 20 iterations. This should NEVER happen. 
  k1=find(k);ix=1;
  while(any(k1) && ix <= count_limit), 
    if trace,figure(ix),end % for debugging 
    [tmpy(k1), tmpt(k1)]= gaussq('tay81fun',uk_inf(k1,ix),uk_inf(k1,ix+1),...
				 tol/2/count_limit,trace,gn,xk(k1)./ak(k1),bk(k1));
    
    yk(k1)=yk(k1)+tmpy(k1);       
    tol1(k(k1))=tol1(k(k1))+tmpt(k1);
     %find integrals which did not converge  
    k1=find(tol1(k) < abs(tmpy) & tol/count_limit <2*abs(tmpy)); 
   
    %k=find(tol1 > abs(tol*int)| tol1 > abs(tol));
    %indices to integrals which did not converge
    ix = ix + 1;
  end
  %yk=yk.*xk./(ak.^2);
  
  k2=find(yk<0 );
  if any(k2)
    yk(k2)= zeros(size(k2));
  end

  k3=find(xk>ak.*sqrt(bk) );
  if any(k3)
    yk(k3)= yk(k3).*(1-4/pi.*acos(ak(k3).*sqrt(bk(k3))./xk(k3))); 
  end
  y(k)=yk.*xk./ak.^2;
end

k4=find(a > 0 & x >= 0 & b==inf);
if any(k4), % special case b==inf -> rayleigh distribution
 y(k4) = pdfray(x(k4),a(k4)/sqrt(2));
end 
y=reshape(y,N,M);


function z=traylfunzeros(x,a,b,n)
% internal function which calculates the approximate zerocrossings of 
%  the integrand:  u* J_o(u/sqrt(b))^b*J_o(x/a* u)
nx=1:n;
z1=a./x; 
z0=sqrt(b);
k=(mod(b,2)==0); % remove
z0(k)=pi/2*(n+1)*max(z1);%size(x)
z0=z0(:,ones(1,n)).* pi.*(nx(ones(size(x)),:)-0.25);
z1=z1(:,ones(1,n)).* pi.*(nx(ones(size(x)),:)-0.25);
z=sort([z0,z1],2);
return
