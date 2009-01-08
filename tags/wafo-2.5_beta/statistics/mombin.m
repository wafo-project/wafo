function [m,v,sk,ku]= mombin(varargin)
%MOMBIN Mean and variance for the BINOMIAL distribution.
% 
% CALL:  [m,v,sk,ku] = mombin(n,p)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   n, p = parameters of the binomial distribution.
%
%  Mean (m) and variance (v) for the Binomial distribution is
%
%  m=n*p  and  v=n*p*(1-p);
%
% Example:
%   par = {10,0.2}
%   X = rndbin(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]        % Estimated mean and variance
%   [mom{1:4}] = mombin(par{:}) % True mean and variance
%
% See also pdfbin, cdfbin, invbin, rndbin, fitbin


Np = 2;
error(nargchk(1,Np,nargin))
options = []; %struct; % default options
params = parsestatsinput(Np,options,varargin{:});

[n,p] = deal(params{:});
if isempty(p)
  error('Probability p undefined!')
end


try

  p(p<0 | p>1) = nan;
  m=n.*p;
  v=n.*p.*(1-p);
  if nargout>2
    sk = (1-2*p)./sqrt(v);
    ku = 3+(1-6.*p.*(1-p))./v;
  end
catch
  error('n and p must be of common size or scalar.');
end