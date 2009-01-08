function f = pdfgengam(x,varargin)
%PDFGENGAM Generalized Gamma probability density function
%
% CALL:  f = pdfgengam(x,a,b,c,options);
%        f = pdfgengam(x,phat,options);
%
%        f = density function evaluated at x
%    a,b,c = parameters   (default b=1,c=1)
%     phat = Distribution parameter struct
%            as returned from FITGENGAM.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
% The generalized Gamma distribution is defined by its pdf
%
% f(x;a,b,c)=c*x^(a*c-1)/b^(a*c)*exp(-(x/b)^c)/gamma(a), x>=0, a,b,c>0
% 
% Example: 
%   x = linspace(0,7,200);
%   p1 = pdfgengam(x,1,2,1); p2 = pdfgengam(x,3,1,1);
%   plot(x,p1,x,p2), shg
%
% See also  cdfgengam, invgengam, rndgengam, fitgengam, momgengam

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History:
% revised pab Aug 2007
% removed call to comnsize-> simpler and faster
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 09.08.2000


options = struct('logp',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  f = options;
  return
end

error(nargchk(2,inf,nargin))
Np = 3;

[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});
if isempty(b),  b=1; end
if isempty(c),  c=1; end

a(a<=0) = nan;
b(b<=0) = nan;
c(c<=0) = nan;
x(x<0)  = inf; % trick to set pdf->0

try 
  xn = (x./b).^c;
  ac = a.*c;
  
  msgID = 'MATLAB:log:logOfZero';
  state  = warning('off',msgID);
  logx =  min(log(x),100*log(realmax));
  logb =  log(b);
  warning(state)
  gammalna = gammaln(a);
  
  %  f=c.*x.^(a.*c-1)./b.^(a.*c)....
  %      .*exp(-(x./b).^c)./gamma(a);
  logfc = (ac-1).*logx-xn-gammalna-ac.*logb;
  %f = c.*exp((ac-1).*logx-xn-gammalna-ac.*logb);
  
  
  % Handle special cases for x == 0 (or logx==-inf):
  spcase = (logx==-inf ) & (ac<=1) & (b>0);
  if any(spcase)
    limit1 = (-gammalna-logb);
    if ~isscalar(limit1)
      limit1 = limit1(spcase & ac==1);
    end
    logfc(spcase & ac==1) = limit1;
    logfc(spcase & ac<1) = inf;
  end  
  if options.logp
    f = logfc+log(c);
  else
    f = c.*exp(logfc);
  end
  
catch
  error('x, a, b and c must be of common size or scalar.');
end
