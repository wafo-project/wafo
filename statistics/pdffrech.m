function f = pdffrech(x,varargin)
%PDFFRECH Frechet probability density function
%
% CALL:  f = pdffrech(x,a,c,options);
%        f = pdffrech(x,phat,options);
%
%        f = density function evaluated at x
%     a, c = parameters
%     phat = Distribution parameter struct
%            as returned from WFRECHFIT.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%
% The Frechet distribution is defined by its cdf
%
%  F(x;a,c) = exp(-(x/a)^(-c)), x>=0, a,c>0
%
% Example: 
%   x = linspace(0,6,200);
%   p1 = pdffrech(x,1,1); p2 = pdffrech(x,2,2); p3 = pdffrech(x,2,5);
%   plot(x,p1,x,p2,x,p3), shg
%
% See also cdffrech, invfrech, rndfrech, fitfrech, momfrech 

% Reference: 

% Tested on; Matlab 5.3
% History: 
% Revised pab Aug 2007
% - removed dependence on comnsize -> calculations faster
% - fixed a bug:  pdffrech(0,a,c) returned nan corrected to infinity for
%                 c<1 and zero otherwise
% Added PJ 10-May-2001


error(nargchk(3,5,nargin))
Np = 2;
options = struct('logp',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a = params{1};
c = params{2};



c(c<=0) = nan;
a(a<=0) = nan; 
x(x<0)  = 0;
try
 

  msgID = 'MATLAB:log:logOfZero';
  state = warning('off',msgID);
  
  logxc = -log(x./c);
  xn    = exp(-c.*log(x./a));
  
  warning(state)
  
  % Make sure f(x==0) = inf if c<1 and f(x==0) = 0 if c>=1
  infty      = realmax;
  islogxcinf = logxc==inf;
  logxc(islogxcinf) = -infty*(c>=1) + infty*(c<1);

  xn(xn==inf & islogxcinf) = 1;
  if options.logp
    f = logxc - xn + log(xn);
  else
    f = xn.*exp(logxc - xn);
  end
  %xn = (a./x).^c;
  %f = xn.*c./x.*exp(-xn);
catch
  error ('x, a and c must be of common size or scalar');
end
return
