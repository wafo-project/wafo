function R = rndinvnorm(varargin)
%RNDINVNORM Random matrices from a Inverse Gaussian distribution.
%
% CALL:  R = rndinvnorm(m0,l,sz);
%        R = rndinvnorm(phat,sz);
%
%  %  R     = a matrix of random numbers from the inverse Gaussian
%            distribution
%      m0,l = parameters  (see pdfinvnorm)
%      phat = Distribution parameter struct
%             as returned from FITINVNORM.  
%        sz = size(R)    (Default common size of m0 and l)
%             sz can be a comma separated list or a vector 
%             giving the size of R (see zeros for options).
%
% Examples:
%   R = rndinvnorm(2,2,100,2);
%   R2 = rndinvnorm(2,3,[100,2]);
%   plotqq(R(:,1),R2(:,1)),shg
%
% See also pdfinvnorm, cdfinvnorm, invinvnorm, fitinvnorm, mominvnorm

% Reference: Chhikara & Folks, "The Inverse Gaussian Distribution", p. 53

% Tested on; Matlab 5.3
% History:
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 14.08.2000

error(nargchk(1,inf,nargin))
Np = 2;
options = struct;
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

m0 = params{1};
l  = params{2};

if isempty(rndsize)
  [iscmn m0 l] = iscomnsize(m0,l);
else
  [iscmn m0 l] =iscomnsize(m0,l,zeros(rndsize{:}));
end
if ~iscmn
  error('m0 and l must be of common size or scalar.');
end
R=zeros(size(m0));
ok=((m0>0)&(l>0));
k=find(ok);
if any(k)
  ksiz=size(k);
  R1=rand(ksiz);
  %Y=rndchi2(1,ksiz);
  Y = rndgam(1/2,2,ksiz);
  X1=m0(k)./(2*l(k)).*(2*l(k)+m0(k).*Y-(4*l(k).*m0(k).*Y+m0(k).^2.*Y.^2).^(1/2));
  X2=m0(k).^2./X1;
  I=(R1<m0(k)./(m0(k)+X1));
  R(k)=X1.*I+X2.*(1-I);
end
k1=find(~ok);
if any(k1)
  R(k1)=nan;
end


