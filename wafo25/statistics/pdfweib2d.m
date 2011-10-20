function f = pdfweib2d(x1,x2,varargin)
%PDFWEIB2D 2D Weibull probability density function (pdf).
%
%  CALL:  f = pdfweib2d(x1,x2,A1,B1,A2,B2,C12,options)
%         f = pdfweib2d(x1,x2,phat,options)
%
%     f    = PDF evaluated at x1 and x2.
%    x1,x2 = coordinates
%   A1, A2 = scale parameters    
%   B1, B2 = shape parameters
%      C12 = interaction parameter between X1 and X2
%     phat = Distribution parameter struct
%            as returned from FITWEIB2D.  
%  options = struct with fieldnames:
%     .logp    : if TRUE, density, p, returned as log(p).
%     .condon  : 0 it returns the the regular pdf of X1 and X2 (default)
%                1 it returns the conditional pdf given X1
%                2 it returns the conditional pdf given X2
%     .meshgrid: if TRUE output computed on meshgrid(x1,x2) (default false)
%     .wdata   : if TRUE output f as wdata object (default false)
%
%   The size of f is the common size of the input arguments X1 and X2. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%   The PDF is defined by:
%
%  f(X1,X2)=B1*B2*xn1^(B1-1)*xn2^(B2-1)/A1/A2/N*...
%           exp{-[xn1^B1 +xn2^B2 ]/N }*I0(2*C12*xn1^(B1/2)*xn2^(B2/2)/N) 
%   where 
%      N=1-C12^2, xn1=X1/A1,  xn2=X2/A2 and 
%      I0 is the modified  bessel function of zeroth order.
%
%  (The marginal distribution of X1 is weibull with parameters A1 and B1) 
%  (The marginal distribution of X2 is weibull with parameters A2 and B2) 
%  (C12 is the interaction parameter between X1 and X2)
%
%Example: 
%  x = linspace(0,6,200); [X1 X2]=meshgrid(x);
%  f = pdfweib2d(X1,X2,2, 2,  3, 2.5, .8);
%  mesh(x,x,f),
%  figure(2), contour(x,x,f)
%
% See also  pdfweib, cdfweib2d, rndweib2d, fitweib2d, momweib2d, besseli

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

 

%  References:
%  Dag Myrhaug & Håvard Rue
%  Journal of Ship Research, Vol 42, No3, Sept 1998, pp 199-205 

%Tested on: matlab 5.2
% History:
% revised pab okt 2007
% - removed dependence on comnsize
%  by Per A. Brodtkorb 13.11.98


Np = 5;
options = struct('logp',false,'condon',0,'meshgrid',false,'wdata',false); % default options
if nargin==1 && strncmpi(x1,'defaults',6)
  f = options;
  return
end
error(nargchk(3,15,nargin))
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a1,b1,a2 b2,c12] = deal(params{:});

if isempty(a1)||isempty(b1)||isempty(a2)||isempty(b2)||isempty(c12)
  error('Requires either 7 input arguments or that input argument 3 is FDATA.'); 
end
if options.meshgrid
  [X1, X2]=meshgrid(x1,x2);
  f1 = localpdfweib2d(X1,X2,a1,b1,a2,b2,c12,options);
else
  f1 = localpdfweib2d(x1,x2,a1,b1,a2,b2,c12,options);
end
if options.wdata
  f = wdata(f1,{x1,x2});
  [cl, pl]=qlevels(f1);
  f = set(f,'contourLevels',cl,'percentLevels',pl,'note','2D Weibull');  
else
  f = f1;
end


function f = localpdfweib2d(x1,x2,a1,b1,a2,b2,c12,options)

a1(a1<=0) = nan;
b1(b1<=0) = nan;
a2(a2<=0) = nan;
b2(b2<=0) = nan;
c12(abs(c12)>=1) = nan;
try
  N1 = 1-c12.^2;
  xn1 = x1./a1;
  xn1(xn1<0) = inf; % trick to set pdf to zero
  xn2 = x2./a2;
  xn2(xn2<0) = inf; % trick to set pdf to zero
  
  zk = 2*c12.*(xn1.^(b1/2).*xn2.^(b2/2))./N1;
  zmax = 32767;
  
  if any(inf>zk(:) & zk(:)>zmax)
    warning('WAFO:PDFWEIB2D','Numerical problems expected!')
  end
  zk = min(zk,zmax);
  [y, ierr] =besseli(0,min(zk,zmax),1 );
  if any(ierr(:)>0)
    y(ierr==1) = nan;
    if any(ierr(:)==2), disp('Overflow.'), end
    if any(ierr(:)==3), disp('Some loss of accuracy in argument reduction.'), end
    y(~isnan(zk) & (ierr>3 | ierr==2)) = 0.0022;
  end
  y = min(log(y),realmax);
catch
   error('Requires non-scalar arguments to match in size.');
end
switch options.condon
  case 0
    limit1 = log(b1./a1);
    logf1 = -xn1.^b1./N1 +log(min(xn1.^(b1-1),realmax)) + limit1;
    
  
    limit2 = log(b2./a2);
    logf2 = -xn2.^b2./N1 +log(min(xn2.^(b2-1),realmax)) + limit2;
    
    logf = logf1+logf2-log(N1)+abs(zk)+y;
    spcase = ( x1==0 & x2 == 0 & b1 < 1& b2+b1-2<0);
    logf(spcase) = inf;
  case 1, %pdf conditioned on x1 ie. p(x2|x1)
    logf = - (c12.*xn1.^(b1/2) - xn2.^(b2/2)).^2./N1 +y  ...
      +log(b2./(a2.*N1))+log(min((xn2).^(b2 - 1),realmax));
    
  case 2,%pdf conditioned on x2 ie. p(x1|x2)
    logf = - (xn1.^(b1/2) - c12.*xn2.^(b2/2)).^2./N1 +y  ...
      +log(b1./(a1.*N1))+log(min((xn1).^(b1 - 1),realmax));
  case 3, % secret option  used by momweib2d: returns x1*p(x1|x2)
    logf = - (xn1.^(b1/2) - c12.*xn2.^(b2/2)).^2./N1 +y  ...
      +log(b1./N1)+b1.*log(min(xn1,realmax));
  case 4, % secret option  used by momweib2d: returns x1^2*p(x1|x2)
     logf = - (xn1.^(b1/2) - c12.*xn2.^(b2/2)).^2./N1 +y  ...
      +log(b1.*a1./N1)+(b1+1).*log(min(xn1,realmax));
  otherwise , error('Illegal value for CONDON')
end

if options.logp
  f = logf;
else
  f = exp(logf);
end
    
return

% Old call kept just in case
% [icode, x1, x2, a1,b1,a2, b2, c12] = iscomnsize(x1,x2,a1,b1,a2,b2,c12);
% if ~icode 
%   error('Requires non-scalar arguments to match in size.');
% end
% xs=size(x1); 
% 
% y = zeros(xs);
% 
% ok = ((a1 > 0) .*(b1 > 0).*(a2 > 0).*(b2 > 0).*(abs(c12)<1));
%   
% k1 = find(~ok);
% if any(k1)
%    tmp   = NaN;
%    y(k1) = tmp(ones(size(k1)));
% end
% 
% k = find((x2 > 0).*(x1 > 0) & ok);
% if any(k),
%   %scales bessel with exp(-abs(x))  to avoid overflow 
%   zk = 2*c12(k).*((x1(k)./a1(k)).^(b1(k)/2).*....
%       (x2(k)./a2(k)).^(b2(k)/2))./(1-c12(k).^2);
%   [y(k), ierr] =besseli(0,zk,1 );
%   
%   switch ierr(1),
%     case 0, %computation OK
%     case 1, error('Illegal arguments.')
%     case 2, disp('Overflow.  Return Inf.')
%     case 3, disp('Some loss of accuracy in argument reduction.')
%     case 4, error('Complete loss of accuracy, z or nu too large.')
%     case 5, error('No convergence.  Return NaN.')
%   end
%   y(k)=log(y(k));
%   
%   switch condon
%     case 0,%regular pdf
%       y(k)=exp(- ((x1(k)./a1(k)).^b1(k) +(x2(k)./a2(k)).^b2(k) ...
% 	  -abs(2*c12(k).*((x1(k)./a1(k)).^(b1(k)/2).*...
% 	  (x2(k)./a2(k)).^(b2(k)/2))))./(1-c12(k).^2)  +y(k) ) ...
% 	  .*b1(k).*(x1(k)./a1(k)).^(b1(k)-1)./a1(k).*b2(k).*...
% 	  (x2(k)./a2(k)).^(b2(k) - 1) ./a2(k)./(1-c12(k).^2);
%     case 1, %pdf conditioned on x1 ie. p(x2|x1)     
%       y(k) =exp(- (c12(k).*(x1(k)./a1(k)).^(b1(k)/2) - ...
% 	   (x2(k)./a2(k)).^(b2(k)/2)).^2./(1-c12(k).^2 ) +y(k)  )...
% 	   .*(b2(k)./(a2(k).*(1-c12(k).^2))).*(x2(k)./a2(k)).^(b2(k) - 1);
%     case 2,%pdf conditioned on x2 ie. p(x1|x2)
%         y(k)=exp(- ((x1(k)./a1(k)).^(b1(k)/2) -c12(k).*...
% 	    (x2(k)./a2(k)).^(b2(k)/2)).^2./(1-c12(k).^2)  +y(k) )...
% 	  .*(b1(k)./(a1(k).*(1-c12(k).^2))).*(x1(k)./a1(k)).^(b1(k)-1);
%     case 3, % secret option  used by momweib2d: returns x1*p(x1|x2) 
%       y(k)=x1(k).*exp(- ((x1(k)./a1(k)).^(b1(k)/2) -c12(k).*...
% 	  (x2(k)./a2(k)).^(b2(k)/2)).^2./(1-c12(k).^2)  +y(k) )...
% 	  .*(b1(k)./(a1(k).*(1-c12(k).^2))).*(x1(k)./a1(k)).^(b1(k)-1);
%     case 4, % secret option  used by momweib2d: returns x1^2*p(x1|x2) 
%       y(k)=x1(k).^2.*exp(- ((x1(k)./a1(k)).^(b1(k)/2) -c12(k).*...
% 	  (x2(k)./a2(k)).^(b2(k)/2)).^2./(1-c12(k).^2)  +y(k) )...
% 	  .*(b1(k)./(a1(k).*(1-c12(k).^2))).*(x1(k)./a1(k)).^(b1(k)-1);
%     otherwise , error('Illegal value for CONDON')
%   end
% 
%   kc=find(isinf(y(k)));
%   if any(kc), 
%     disp('Computational problem occured: returning NaNs') ;   
%     y(k(kc))=NaN;
%   end
% end
% 
% % Special case for asymptotes.
% switch condon
%   case 0,k1 = find( (x1 == 0 & b1 < 1& x2>0)|( x2 == 0 & b2 < 1& x1>0 )|  ...
% 	( x1==0 & x2 == 0 & b1 < 1& b2+b1-2<0)    );
%   case 1, k1 = find( x2 == 0 & b2 < 1 );
%   case {2},  k1 = find( x1 == 0 & b1 < 1 );
%   case {3,4}, k1=[]; %do nothing
% end
% if any(k1)
%   tmp   = Inf;
%   y(k1) = tmp(ones(size(k1)));
% end
% 
% % Special case when the marginal Weibull is the same as exponential. 
% switch condon,
%   case 0, 
%     k2 = find((x1 == 0 & b1 == 1&x2>0)  );
%     if any(k2),
%       y(k2) =  b2(k2) .* (x2(k2)./a2(k2)).^ (b2(k2) - ...
% 	  1)./a2(k2)./(1-c12(k2).^2)./a1(k2).*....
% 	  exp(- (  (x2(k2)./a2(k2)) .^ b2(k2) )./(1-c12(k2).^2)   );
%     end
%     k3 = find((x2 == 0 & b2 == 1&x1>0)  );
%     if any(k3),
%       y(k3) = b1(k3) .* (x1(k3)./a1(k3)).^ (b1(k3) - ...
% 	  1)./a1(k3)./a2(k3)./(1-c12(k3).^2).*...
% 	  exp(- ((x1(k3)./a1(k3)) .^ b1(k3)  )./(1-c12(k3).^2)   );
%     end
%     k4 = find( x1==0 & x2==0 & b1 == 1& b2==1)   ;
%     if any(k4),
%       y(k4) = 1./(a1(k4).*a2(k4).*(1-c12(k4).^2));
%     end
%     
%   case 1, 
%     k3 = find(x2 == 0 & b2 == 1  ); 
%     if any(k3),
%       y(k3) = exp(- (c12(k3).^2.*(x1(k3)./a1(k3)) .^ b1(k3) ...
% 	  )./(1-c12(k3).^2)   )./a2(k3)./(1-c12(k3).^2) ;
%     end
%     
%   case {2}, 
%     k3 = find(x1 == 0 & b1 == 1  ); 
%     if any(k3),
%       y(k3) = exp(- (c12(k3).^2.*(x2(k3)./a2(k3)) .^ b2(k3) ...
% 	  )./(1-c12(k3).^2)   )./a1(k3)./(1-c12(k3).^2) ;
%     end
%   case {3,4}, %do nothing
% end 