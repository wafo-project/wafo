function [r]=ratio(a,b,sa,sb)
%RATIO compute ratio of hyperbolic functions
%      to allow extreme variations of arguments.
%
% CALL: r=ratio(a,b,sa,sb);
%
%       r = f(a,sa)./f(b,sb), ratio vector hyperbolic functions of same
%                size as a and b
%     a,b = arguments vectors of the same size
%   sa,sb = defines the hyperbolic function used, i.e.,
%           f(x,1)=cosh(x), f(x,-1)=sinh(x)
%
% Examples:
% ratio(2,1,1,1)   % gives r=cosh(2)/cosh(1)
% ratio(2,1,1,-1)  % gives r=cosh(2)/sinh(1)
% ratio(2,1,-1,1)  % gives r=sinh(2)/cosh(1)
% ratio(2,1,-1,-1) % gives r=sinh(2)/sinh(1)
%
% See also  tran

% Tested on: matlab 5.2
% history
% revised pab dec2003 
% commented out old call  
% added todo comment  
% revised pab 09.10.2002
% -fixed bug: replaced * with .* thanks to Françoise GIRARD 
% -added more checks when a==b and when a<0 or b<0 => made it more robust
% revised pab 07.11.2001
% -added comnsize + see also line
% -Fixed a bug: ratio(0,0,-1,-1) gave NaN but should return 1
% revised pab 11.01.2000
% - added sign(s1), sign(s2) to ensure correct calculation
% - updated documentation
% - fixed a bug in expression.
% by L. Borgman

% TODO % Does not always handle division by zero correctly

[iscmn, a,b,sa,sb]=iscomnsize(a,b,sign(sa),sign(sb));
if ~iscmn,
   error('a,b,sa and sb must be of common size or scalar!')
end	
r = zeros(size(a));

sc = (a==b);
 
k = find(sc);
if any(k)
   r(k) = 1;
   d = 0.5*(sb(k)-sa(k));
   k0 = find(d~=0);
   if any(k0)
     k00 = k(k0); 
     r(k00) = tanh(a(k00)).^d(k0);
   end	
end	


k1 = find(~sc);
if any(k1),
   ak = a(k1);
   bk = b(k1);
   sak = sa(k1);
   sbk = sb(k1);
   signRatio = ones(size(k1));
   ka = find(sak.*ak<0);
   if any(ka)
     signRatio(ka) =  sak(ka);
  end	
  kb = find(sbk.*bk<0);
  if any(kb)
     signRatio(kb) =  signRatio(kb).*sbk(kb);
  end	

  %signRatio = (2*(0<=ak)-1).*(2*(0<=bk)-1).*sak.*sbk
   bk = abs(bk);
   ak = abs(ak);
   msgId = 'MATLAB:divideByZero';
   state = warning('off',msgId);   % fix to avoid warning messages about division by zero.    
   r(k1)=signRatio.*exp(ak-bk).*((sak+exp(-2*ak))./(sbk+exp(-2*bk)));
   warning(state)
 end	

return

% Old call
%r=exp(a-b).*(1+sign(sa)*exp(-2*a))./(1+sign(sb)*exp(-2*b));


 
