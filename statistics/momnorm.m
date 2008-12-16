function [m,v,sk,ku]= momnorm(varargin)
%MOMNORM Mean and variance for the Normal distribution.
% 
% CALL:  [m,v] = momnorm(m0,v0)
%
%   m, v = the mean and variance, respectively 
% m0, v0 = parameters of the Normal distribution.
%
%  Mean (m) and variance (v) for the Normal distribution is
%
%  m=m0  and  v=v0;
%
% Example:
%   par = {-1,1}
%   X = rndnorm(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]       % Estimated mean and variance
%   [m,v] = momnorm(par{:}) % True mean and variance
%
% See also pdfnorm, cdfnorm, invnorm, rndnorm, fitnorm
Np = 2;
error(nargchk(1,Np,nargin))
options = [];
params = parsestatsinput(Np,options,varargin{:});

[m0,v0] = deal(params{:});
if isempty(v0)
  v0 = 1;
end
if isempty(m0)
  m0=0;
end
v0(v0<=0) = nan;
m =  m0;
v = v0;
sk = 0;
ku = 3;

[iscmn,m,v,sk,ku ] = iscomnsize(m,v,sk,ku);
if ~iscmn
  error ('m and v must be of common size or scalar');
end
