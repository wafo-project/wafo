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
%   par = {10,0.2};
%   X = rndbin(par{:},10000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = mombin(par{:}); % True mean and variance
%   assert(moments, mom, 0.05)
%
% See also pdfbin, cdfbin, invbin, rndbin, fitbin

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



Np = 2;
%error(nargchk(1,Np,nargin))
narginchk(1,Np)
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