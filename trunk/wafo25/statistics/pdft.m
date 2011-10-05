function f = pdft(x,varargin)
%PDFT  Student's T probability density function
%
% CALL:  f = pdft(x,df)
%
%   f = density function evaluated at x
%   x = matrix
%   df = degrees of freedom (1,2,....)
%
% The Cauchy distributions equals the T distribution when df=1.
% As df->inf the T distribution converges to the normal distribution.
%
% Example:
%   x = linspace(-5,5,200);
%   p1 = pdft(x,1); p2 = pdft(x,5);
%   plot(x,p1,x,p2,x,pdfnorm(x)),shg
%
% See also cdft, invt, rndt, fitt, momt

% tested on matlab 5.3
%History:
% revised pab oct 2007
% -new parametrization which is more stable for large df.
%revised pab 29.10.2000
% adapted from stixbox changed name to pdft
% -added nargchk + check on floor(df)==df
% - changed from gamma to gammaln for more stable computation
% - added the secret option disable in order to use this function for MLE
%   estimation 
%  by Anders Holtsberg, 18-11-93
%     Copyright (c) Anders Holtsberg

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

options = struct('logp',false,'disable',false); % default options
if (nargin==1 && nargout <= 1 && isequal(x,'defaults'))
  f = options; 
  return
end

error(nargchk(2,inf,nargin))
Np = 1;

[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
mxdf = 1e20;
df = params{1};
df(df>mxdf) = inf;
%df = min(params{1},1e20);
%df = min(params{1},realmax);

try 
  if options.disable 
    df(df<=0|df~=floor(df)) = nan;
  else
    df(df<=0) = nan;
  end
  %mxdf = 10^8;
  % const = gammaln((df+1)/2)-gammaln(df/2)-0.5*log(df/2) ;
  % const = gammaln((df+3)/2)-gammaln((df+2)/2)+log(sqrt(2*df)/(df+1));
  % const(df>mxdf) = 0;
  const = (stirlerr((df+1)/2)-stirlerr(df/2)+(df.*log1p(1./df)-1)/2);
  const(df>mxdf) = 0;
  x2 = x.^2;
  corr1 = (1+1./df).*log1pxdx(x2./df);
  
  if options.logp
    f = const-0.5*x2.*corr1-0.5*log(2*pi);
  else
    f = exp(const-0.5*x2.*corr1)./sqrt(2*pi);
  end
catch
   error('x and df must be of common size or scalar');
end



