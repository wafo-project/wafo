function [F,Flo,Fup] = cdft(x,varargin)
%CDFT  Student's T  cumulative distribution function
%
% CALL:  F = cdft(x,df,options);
%        [F,Flo,Fup] = cdft(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        x = matrix
%       df = degrees of freedom (1,2,....)
%     phat = Distribution parameter struct
%            as returned from FITT.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% Example:
%   x = linspace(-5,5,200);
%   p1 = cdft(x,1); p2 = cdft(x,5);
%   plot(x,p1,x,p2), shg
%
% See also pdft, invt, rndt, fitt, momt

% tested on matlab 5.3
%History:
%revised pab 22.05.2003
% -added new methods for df==1 or df==2 and for region1= x^2<df^2
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize, mxdf +  check on floor(df)==df
%by      Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg


options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'disable',false); % default options
if (nargin==1 && nargout <= 1 && isequal(x,'defaults'))
  F = options; 
  return
end

error(nargchk(2,inf,nargin))


Np = 1;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
df = params{1};

[icode x,df] = iscomnsize(x,df);
if ~icode
  error('x and df must be of common size or scalar');
end

F = zeros(size(x));

%df = min(df,1000000); % make it converge and also accept Inf.
mxdf = 10^6;
ok      = (0<df);

if ~options.disable || all(df==floor(df))
  ok = ok & df==floor(df);
  k1 = find(df==1 );
  if any(k1)
    F(k1) =  ( 1 + 2*atan(x(k1))/pi )/2;
  end
  k2 = find(df==2);
  if any(k2)
    x2= x(k2);
    F(k2) = ( 1 + x2./sqrt( 2 + x2.*x2 ))/2;
  end
  
  region0 = (2<df) & df<mxdf ;
  region1 = (abs(x)<sqrt(abs(df))) ;
  region2 = region0 & ~region1;
  k3 = find(ok & region1 & region0);
  if (any(k3)),
    dfk = df(k3);
    xk = x(k3);
    xk2 = xk.*xk;
    cssthe = 1./( 1 + xk2./dfk );
    nuVec = unique(dfk(:));
    Fk = zeros(size(dfk));
    for nu = nuVec(:).'
      knu = find(dfk==nu);
      Fk(knu) = evalPoly(nu,xk(knu),xk2(knu),cssthe(knu));
    end
    F(k3) = Fk;
  end
else 
  region2 = df<mxdf;
end

k = find(ok & region2);
if any(k),
%   sgn = sign(x(k)) + x(k)==0;
%   
%   probF = cdff(x(k).^2,1,df(k),'disable',options.disable);
%  
%   F(k) = (1 + sgn.*probF)/2;h
  neg = x(k)<0;
  %tmp = 1-(1-cdff(x(k).^2,1,df(k)))/2;
  tmp = (1+cdff(x(k).^2,1,df(k),'disable',options.disable))/2;
  F(k) = tmp + (1-2*tmp).*neg;
  %xn = 1./(1+df(k)./(x(k).^2));
  %tmp = cdfbeta(xn,1/2,df(k)/2,'lowertail',false)/2;

end

k1=find(ok & df>=mxdf);
if any(k1)
  F(k1) = cdfnorm(x(k1),0,1);
end

  
k2 = find(~ok);
if any(k2)
  F(k2)=NaN;
end

if nargout>1
% TODO % Implement  Flo and Fup
 warning('WAFO:PRBT','Flo and Fup not implemented yet')
 Flo = nan;
 Fup = Flo;
end



if options.logp
  if options.lowertail
    F = log(F);
  else
    F = log1p(-F);
  end
elseif ~options.lowertail
 F = 1-F;
end
  

return
function F = evalPoly(nu,t,tt,cssthe)
  
  polyn = 1;
  for j = nu-2 : -2 : 2
    polyn = 1 + ( j - 1 )*cssthe.*polyn/j;
  end 
  if ( mod( nu, 2 ) == 1 ) 
    ts = t/sqrt(nu);
    F = ( 1 + 2*( atan(ts) + ts.*cssthe.*polyn )/pi )/2;
  else
    snthe = t./sqrt( nu + tt );
    F = ( 1 + snthe.*polyn )/2;
  end
  F = max(0, min(F, 1) );
  return