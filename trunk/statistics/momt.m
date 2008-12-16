function [m,v,sk,ku]= momt(df)
%MOMT Mean and variance for the Student's T  distribution.
% 
% CALL:  [m,v,sk,ku] = momt(df)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   df   = degrees of freedom of the Student's T distribution
%
%  Mean (m) and variance (v) for the T distribution is
%
%  m=0 if df>1  and  v=df/(df-2) if df>2
%
% Example:
%   par = {10}
%   X = rndt(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]       % Estimated mean and variance
%   [m,v,sk,ku] = momt(par{:}) % True mean and variance
%
% See also pdft, cdft, invt, rndt, fitt



% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% by pab 23.10.2000

error(nargchk(1,1,nargin))


m = zeros(size(df));
m(df<=1) = nan;
df(df<=2) = nan;
v = df./(df-2);
df(df<=3) = nan;
sk = m;
sk(isnan(df)) = nan;
df(df<=4) = nan;
ku = 3+6./(df-4);


