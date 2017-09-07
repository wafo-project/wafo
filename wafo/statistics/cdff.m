function [F,Flo,Fup] = cdff(x,varargin)
%CDFF  Snedecor's F cumulative distribution function
%
% CALL:  F = cdff(x,df1,df2,options);
%        [F,Flo,Fup] = cdff(x,df1,df2,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        x = matrix
%  df1,df2 = degrees of freedom (1,2,....)
%     phat = Distribution parameter struct
%            as returned from WFFIT.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .disable  : if TRUE disable check on integer DF.
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The F- CDF is defined by
% 
%   F(x;df1,df2) = betainc(1/(1+df2/(df1*x),df1/2,df2/2)
%
% Example:
%   x  = linspace(0,6,200);
%   p1 = cdff(x,1,1); p2 = cdff(x,2,2);
%   plot(x,p1,x,p2), shg
%
% See also pdff, invf, rndf, wffit, momf

% tested on matlab 5.3
%History:
% -revised pab oct 2007
% - more robust for large df1 and df2
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize +  check on floor(df)==df
%        Anders Holtsberg, 18-11-93
%        Copyright (c) Anders Holtsberg

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


options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'disable',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
%error(nargchk(2,9,nargin))
narginchk(2,9)
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

df1 = params{1};
df2 = params{2};

ok1= ((df1 >0) );
ok2= ((df2 >0));
if options.disable
  warntxt =  'df%d should be positive.';
else
   ok1 = (ok1 & df1==round(df1));
   ok2 = (ok2 & df2==round(df2) );
   warntxt =  'df%d should be a positive integer'; 
end

if any(~ok1)
  warning('WAFO:PRBF',warntxt,1)
  df1(~ok1) = nan;
end

if any(~ok2)
  warning('WAFO:PRBF',warntxt,2)
  df2(~ok2) = nan;
end

x(x<0)  = 0;
try
   c = (df2./df1);
  xn  = sign(x).*exp(-log1p(c./abs(x)));
  %xn = 1./(1+df2./(df1.*x));
catch
  error('x, df1 and df2 must be of common size or scalar');
end
if any(abs(df1(:)*eps)>abs(df2(:)) | abs(df1(:))>realmax/2) %any(abs(k(:)*eps)>2) 
  warning('WAFO:PRBF','Result may not be exact.')
end

lrg = sqrt(eps)*c>1 & (x>=0) ;
if nargout>1
  opt = options;
  opt.covariance = options.covariance/4;
  [F,Flo,Fup] = cdfbeta(xn,df1/2,df2/2,opt);
  if any(lrg)
    % c-> inf => a*X -> pdfchi2(a)
    if ~isscalar(x), x = x(lrg);end
    if ~isscalar(df1), df1 = df1(lrg);end
    opt.covariance = options.covariance(1);
    [F(lrg),Flo(lrg),Fup(lrg)] = cdfchi2(df1.*x,df1,options);
  end
else
  F = cdfbeta(xn,df1/2,df2/2,options);
 
  if any(lrg)
    % c-> inf => a*X -> pdfchi2(a)
    if ~isscalar(x), x = x(lrg);end
    if ~isscalar(df1), df1 = df1(lrg);end
    F(lrg) = cdfchi2(df1.*x,df1,options);
  end
end
% 
% [errorcode x,a,b] = comnsize(x,a,b);
% if errorcode>0,
%   error('x, df1 and df2 must be of common size or scalar');
% end
% 
% F = zeros(size(x));
% 
% ok = (a>0 & b>0 & floor(a)==a & floor(b)==b);
% k=find(ok & x>=0 & x<inf);
% if any(k),
%   F(k) = cdfbeta(x(k)./(x(k)+b(k)./a(k)),a(k)/2,b(k)/2);
% end
% 
% k1=find(ok &  x==inf);
% if any(k1),
%   F(k1) =ones(size(k1));
% end
% 
%   
% k2 = find(~ok);
% if any(k2)
%   warning('df1 and df1 must be positive integers.')
%   f(k2)=NaN;
% end
% 
