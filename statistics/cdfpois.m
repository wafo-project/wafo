%CDFPOIS Poisson cumulative distribution function
%
%  CALL:  F = cdfpois(x,L,options)
%         [F,Flo,Fup] = cdfpois(x,phat,options)
% 
%        F = probability of observing x givel L
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        x = 
%        L = mean and variance of the poisson distribuion
%     phat = Distribution parameter struct
%            as returned from FITPOIS.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Poisson probability mass function is defined by:
%
%        f(x) = L^x exp(-L)/x!, 0<=L, x=0,1,2,....
%
%  Example
%  x = 1:10; L = 4;
%  F = cdfpois(x,L);
%  F2 = 1-cdfgam(L,x+1); % is the same
%
% See also pdfpois, invpois, rndpois, fitpois, mompois

function [F,Flo,Fup] = cdfpois(x, varargin)

options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'disable',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
error(nargchk(2,inf,nargin))
Np = 1;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[lb] = deal(params{:});
lb(lb<0) = nan;
x(x<0) = -1;
if ~options.disable
  %x = round(x); % force int
  x(x~=floor(x)) = -1;
end

if options.lowertail
  tail = 'upper';
else
  tail = 'lower';
end
try
  if options.logp
    F = gammaincln(lb,floor(x)+1,tail);
  else
    F = gammainc(lb,floor(x)+1,tail);
  end
catch
  error ('x and L must be of common size or scalar');
end

% TODO % Implement Flo and Fup
Flo = nan;
Fup = nan;
