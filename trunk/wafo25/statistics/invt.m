function [x,xlo,xup] = invt(F,varargin)
%INVT Inverse of the Student's T distribution function
%
% CALL:  x = invt(F,df);
%
%   x   = inverse cdf for Student's T distribution evaluated at F.
%   df  = degrees of freedom
%   
%  Example:
%    F = linspace(0,1,100);
%    x = invt(F,1);
%    plot(F,x)
%
% See also pdft, cdft, rndt, fitt, momt


% tested on matlab 5.3
%History:
%revised pab 22.05.2003
% - added new method for df==2
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize
% - added check on F, df + some code from plotnorm when df>mxdf
%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg

options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false);
if (nargin==1 && nargout <= 1 && isequal(F,'defaults'))
  x = options; 
  return
end
error(nargchk(2,inf,nargin))
Np = 1;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
df = params{1};
xlo = nan;
xup = nan;
% TODO % Implement xlo and xup output
% TODO % Implement better logp and lowertail options
if options.logp
  if options.lowertail
    F = exp(F);
  else
    F = -expm1(F);
  end
elseif ~options.lowertail
  F = 1-F;
end

[iscmn F,df] = iscomnsize(F,df);
if ~iscmn,
  error('F and df must be of common size or scalar');
end
x = zeros(size(F));

mxdf = 10^2;

ok = (df>0 & floor(df)==df & 0 <=F & F<=1  );
region1 = ((0 < F) & (F < 1));

k1 = find((df == 1) & region1 & ok);
if any(k1)
  x(k1) = tan(pi * (F(k1) - 0.5));
end

k2 = find((df == 2)& region1 & ok);
if (any(k2))
  R = (2*F(k2)-1);
  x(k2) = sqrt(2)*R./sqrt(1-R.*R);
end

k = find(ok & region1 &  2 < df  & df<=mxdf);
if any(k),
  s = (F(k)<0.5); 
  tmp = F(k) + (1-2*F(k)).*s;
  tmp2 = invbeta(1-(2*(1-tmp)),1/2,df(k)/2);
  tmp3 = tmp2.*df(k)./((1-tmp2));
  x(k) = (1-2*s).*sqrt(tmp3);
end
k=find(ok & region1 & df>mxdf);
if any(k), % pab 01.11.2000, added from plotnorm 
  x(k) = invnorm(F(k));
  k0=find(df(k)<inf);
  if any(k0)
    k1=k(k0);
    Y=x(k1);
    g1=1/4*(Y.^3+Y);
    g2=1/96*(5*Y.^5+16*Y.^3+3*Y);
    g3=1/384*(3*Y.^7+19*Y.^5+17*Y.^3-15*Y);
    g4=1/92160*(79*Y.^9+776*Y.^7+1482*Y.^5-1920*Y.^3-945*Y);
    x(k1)=Y+g1./df(k1)+g2./df(k1).^2+g3./df(k1).^3+g4./df(k1).^4;
  end
end
 

k2=find(ok&F==0);
if any(k2)
  x(k2)=-inf;
end

k3=find(ok&F==1);
if any(k3)
  x(k3)=inf;
end




k4=find(~ok);
if any(k4)
  x(k4)=NaN;
end

