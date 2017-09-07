function [m,v,sk,ku]= mombeta(varargin)
%MOMBETA Mean and variance for the Beta distribution.
% 
% CALL:  [m,v,sk,ku] = mombeta(a,b)
%
%   m,  v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   a,  b = parameters of the Beta distribution
%
%  Mean (m) and variance (v) for the Beta distribution is
%
%      m = a/(a+b) and v = a*b/(a+b)^2/(a+b+1)  if a>0, b>0
%
% Example
%   par = {1,1};
%   X = rndbeta(par{:},10000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};         % Estimated mean and variance
%   [mom{1:4}] = mombeta(par{:}); % True mean and variance
%
%   assert(moments, mom, 0.05)
%
% See also  pdfbeta, cdfbeta, invbeta, rndbeta, fitbeta

% Copyright (C) 2000, 2007 Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
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


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% revised pab Aug 2007
% -Simplified code removed dependence on comnsize.
% by pab 23.10.2000

Np = 2;
%error(nargchk(1,Np,nargin))
narginchk(1,Np)
options = [];
params = parsestatsinput(Np,options,varargin{:});

[a,c] = deal(params{:});
if isempty(c)
  error('Shape parameter b undefined!')
end


try
  a(a<=0) = nan;
  c(c<=0) = nan;
  cdap1 = (c+1)./(a+1);
  cda = c./a;
  apc = a+c;
  cda(c==a) = 1;
  cdap1(c==a) = 1;
  m = 1./(1+cda);
  v = m./(1+1./cda)./(apc+1);  
  sk = 	(2*(c-a))./(2+apc).*sqrt((1+apc)./(a.*c));
  sk(c==a) = 0;
  %ku = 3+(6.*(a.^3+a.^2.*(1-2.*c)+c.^2.*(1+c)-2.*a.*c.*(2+c)))./(a.*c.*(2+a+c).*(3+a+c))
  ku = 3+6.*(1./(cda.*(cdap1+1))+cda.*cdap1./(cdap1+1)-2)./(apc+3);
catch
  error('a and b must be of common size or scalar.');
end

% Saved Old call just in case

% [errorcode, a, c] = comnsize(a,c);
% if errorcode > 0
%     error('a and b must be of common size or scalar.');
% end
% 
% 
% %   Initialize Mean and Variance to zero.
% m = zeros(size(a));
% v = zeros(size(a));
% 
% ok = (a > 0  & c > 0 );
% k = find(ok);
% if any(k)
%   m(k) = a(k)./(a(k)+c(k));
%   v(k) = m(k).*c(k)./(a(k)+c(k))./(a(k)+c(k)+1);
% end
% 
% k1 = find(~ok);
% if any(k1)
%   warning('a and b should be positive')
%   tmp = NaN;
%   v(k1) = tmp(ones(size(k1)));
%   m(k1) = v(k1);
% end
% 
