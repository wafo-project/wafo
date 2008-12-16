function [m,v,sk,ku]= momf(varargin)
%MOMF Mean and variance for the Snedecor's F distribution.
% 
% CALL:  [m,v,sk,ku] = momf(df1,df2)
%
%   m,  v   = the mean and variance, respectively 
%     sk,ku = the skewness and kurtosis, respectively. 
%  df1, df2 = degrees of freedom of the F distribution
%
%  Mean (m) and variance (v) for the F distribution is
%
%      m = df2/(df2-2)                  if df2>2  
% and  
%      v = 2*m^2*(df1+df2-2)/(df2-4)/df1  if df2>4
%
% Example
%   par = {9,9}
%   X = rndf(par{:},10000,1);
%   [mean(X) var(X),skew(X),kurt(X)]         % Estimated mean and variance
%   [mom{1:4}] = momf(par{:}) % True mean and variance
%
% See also pdff, cdff, invf, rndf, fitf


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% by pab 23.10.2000

error(nargchk(1,2,nargin))
Np = 2;
options = [];
params = parsestatsinput(Np,options,varargin{:});
[a,c] = deal(params{:});
if isempty(c)
  error('df2 undefined!')
end
c(c<=2 | floor(c)~=c ) = nan;
a(a<=0 | floor(a)~=a ) = nan;
m = c./(c-2);
try
  c(c<=4) = nan;
  v = 2*m.^2.*(c+a-2)./(c-4)./a;

  c(c<=6) = nan;

  sk = 	(2*(2.*a+c-2))./(c-6).*sqrt((2.*(c-4))./(a.*(c+a-2)));

  c(c<=8) = nan;

  ku =   (12.*(20*c-8.*c.^2+c.^3+44.*a-32.*a.*c+5.*c.^2.*a-22.*a.^2+5.*c.*a.^2-16))./(a.*(c-6).*(c-8).*(a+c-2));
  if isscalar(m)
    m = m(ones(size(a)));
  end
catch
  
%[iscmn m,a] = iscomnsize(m,a);
%if ~iscmn
  error('df1 and df2 must be of common size or scalar.');
end
