function f = pdfchi2(x,varargin)
%PDFCHI2 Chi squared probability density function
%
% CALL:  f = pdfchi2(x,df,options);
%        f = pdfchi2(x,phat,options);
%
%        f = density function evaluated at x
%       df = degrees of freedom
%     phat = Distribution parameter struct
%            as returned from FITCHI2.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%         .disable: if TRUE disable check on integer DF. 
% 
% The Chi squared distribution is defined by its pdf
%   f(x)=x^(df/2-1)*exp(-x/2)/gamma(df/2)/2^(df/2), x>=0, df=1,2,3,...
% The CHI^2 is a special case of the gamma distribution, i.e.:
%   pdfchi2(x,df)=pdfgam(x,df/2,2)
%
% Example: 
% x = linspace(0,7,200);
% p1 = pdfchi2(x,2); p2 = pdfchi2(x,3);
% plot(x,p1,x,p2),shg
%
% See also pdfgam, cdfchi2, invchi2, rndchi2, fitchi2, momchi2

% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 415 ff
% Wiley


% Tested on: Matlab 5.3
% History:
% revised pab aug 2007
% -removed call 2 comnsize
% revised pab 25.10.2000
%  - added comnsize, nargchk
%  - replaced code with a call to pdfgengam -> made maintanence easier
% added ms 15.06.2000

options = struct('logp',false,'disable',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  f = options;
  return
end
%error(nargchk(2,inf,nargin))
narginchk(2,inf)
Np = 1;

[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

p = params{1};

try
  % secret option in order to make ML estimation work
  if options.disable,
    p(p<=0) = nan;
  else
    p(p<=0 | p~=round(p)) = nan;
  end
  b = 2;
  c = 1;
  f = pdfgengam(x,p/2,b,c,options);
catch
  error('x and p must be of common size or scalar.');
end
return
% [errorcode,x,p,b,c] = comnsize(x,p,2,1);
% if errorcode > 0
%   error('x and p must be of common size or scalar.');
% end
% 
% f=zeros(size(x));
% if disable,
%   ok = (p>0);
% else
%   ok = (p==round(p) & p>0);
% end
% k = find(ok);
% if any(k),
%   f(k) = pdfgengam(x(k),p(k)/2,b(k),c(k));
%   %f=x.^(p/2-1).*exp(-x/2)/gamma(p/2)/2^(p/2).*(x>=0);
% end
% 
% k1=find(~ok);
% if any(k1),
%   warning('p should be a positive integer')
%   f(k1)=NaN;
% end






