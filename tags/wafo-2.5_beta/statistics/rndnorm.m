function R = rndnorm(varargin)
%RNDNORM Random matrices from a Normal distribution.
%
% CALL:  R = rndnorm(mu,v,sz);
%
%        mu = mean       (Default 0)
%         v = variance   (Default 1)
%        sz = size(R)    (Default common size of mu and v)
%             sz can be a comma separated list or a vector 
%             giving the size of R (see zeros for options). 
%
% Examples:
%   R  = rndnorm(1,2,100,2,2);
%   R2 = rndnorm(1,2,[100,3,2]);
%   plotnorm([R(:,:,1) R2(:,:,2)])
%
% See also  pdfnorm, cdfnorm, invnorm, fitnorm, momnorm

% tested on: matlab 5.3
% History:
% revised pab 23.10.2000
%  - added default mu,v
%  - added comnsize, nargchk
%  - added greater flexibility on the sizing of R
% by ??

error(nargchk(0,inf,nargin))
options = struct; % default options
Np = 2;
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[mu,v] = deal(params{1:Np});
if isempty(mu),mu=0;end
if isempty(v), v=1;end



if isempty(rndsize)
  csize = comnsize(mu,v);
else
  csize = comnsize(mu,v,zeros(rndsize{:}));
end
if any(isnan(csize))
    error('mu and v must be of common size or scalar.');
end
v(v<0)=nan;
R = randn(csize).*sqrt(v)+mu;
