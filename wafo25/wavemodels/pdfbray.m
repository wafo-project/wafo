function f = pdfbray(x,varargin)
%PDFBRAY Beta Rayleigh PDF of wave heigth
%                     
%   f = 2*(a+b-1)!/((a-1)! * (b-1)!)*h^(2a-1)*(1-(h/c)^2)^(b-1)/c^(2a)
%                       
%  CALL:   f = pdfbray(h,a,b,c); 
%
%       f = pdf
%       h = waveheigth (0 <= h <= c)
%       a = abs(k1*(k2-k1)/(k1^2-k2)) 
%       b = abs((1-k1)*(k2-k1)/(k1^2-k2)) 
%       c = Hb, breaking wave height approximated by water depth, d.
% where
%      k1 = E(H^2)/Hb^2
%      k2 = E(H^4)/Hb^4
%  E(H^2) = .5*exp(0.00272*(d/g*Tp^2)^(-0.834))*Hm0^2
%  E(H^2) = .5*exp(0.00046*(d/g*Tp^2)^(-1.208))*Hm0^2
%     Hm0 = significant waveheight
%     Tp  = modal period of wave spectrum
%
%    The size of F is the common size of H, A, B and C.  A scalar input   
%    functions as a constant matrix of the same size as the other input.
%
% Example: % Compare with Rayleigh distribution
%  Hm0 = 7;Tp = 11;d = 50; g = gravity;
%  k1  = .5*exp(0.00272*(d/g*Tp^2)^(-0.834))*Hm0^2/d^2;
%  k2  = .5*exp(0.00046*(d/g*Tp^2)^(-1.208))*Hm0^2/d^4;
%  a   = abs(k1*(k2-k1)/(k1^2-k2)); 
%  b   = abs((1-k1)*(k2-k1)/(k1^2-k2));
%  h   = linspace(0,2*Hm0)';
%  plot(h,pdfbray(h,a,b,d),'r',h,pdfray(h,Hm0/2))
%
% See also  pdfbeta


% 
%   Reference:
%       Michel K. Ochi (1998),
%      "OCEAN WAVES, The stochastic approach",
%       OCEAN TECHNOLOGY series 6, Cambridge, pp 279. (pd of peaks to trough) 

%tested on: matlab 5.2
%History:
% revised pab 31.03.2001
% added example
% revised pab 14.10.1999
% updated help header
% by  Per A. Brodtkorb 21.02.99

error(nargchk(2,inf,nargin))
Np = 3;
options = struct('logp',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});
if isempty(c)
  error(nargchk(4,inf,2))
end
if isempty(b)
  error(nargchk(4,inf,1))
end

a(a<=0) = nan;
b(b<=0) = nan;
c(c<=0) = nan;
x(x<0)  = 0;
try
  xn = (x./c).^2;
  %xn(xn>1) = 1;
  if options.logp
    f = pdfbeta(xn,a,b,options) + log(2.*x./c.^2);
 else
   f = 2*pdfbeta(xn,a,b).*x./c.^2;
   end
catch
    error('h, a, b and c must be of common size or scalar.');
end

% return
% [csize, x, a ,b, c] = comnsize(x,a,b,c);
% 
% if any(isnan(csize))
%     error('h, a, b and c must be of common size or scalar.');
% end
% 
% 
% % Initialize Y to zero.
% y=zeros(size(x));
% 
% 
% k=find(a > 0 & x >=0 & b>0 & c>0 & x<=c & ~((x == 0 & a < .5) | (x == c & b < 1)));
% if any(k),
%   xk = x(k); ak = a(k); bk = b(k);ck=c(k);
%  switch 1, % choose between different implementations
%   case 1, 
%          y(k)=2*(xk./ck).^(2*ak-1)./ck.*(1-(xk./ck).^2).^(bk-1).*exp(-betaln(ak,bk));
%  case 2,
%         tmp(k) = (2*ak - 1).*log(xk./ck)-log(ck) + (bk - 1).*log((1 - (xk./ck).^2)) - betaln(ak,bk);
%         y(k) = 2*exp(tmp(k));
%  case 3,
%       y(k)=2*pdfbeta((xk./ck).^2,ak,bk).*xk./ck.^2;
%  end
% end
% 
% % Return Inf for x = 0 and a < 1 or x = 1 and b < 1.
% % Required for non-IEEE machines.
% k2 = find((x == 0 & a < .5) | (x == c & b < 1));
% if any(k2)
%     tmp = Inf;
%     y(k2) = tmp(ones(size(k2))); 
% end
% 
% % Return NaN if A,B or C  is not positive.
% k1 = find(a <= 0| b<=0|c<=0);
% if any(k1) 
%     tmp   = NaN;
%     y(k1) = tmp(ones(size(k1)));
% end
