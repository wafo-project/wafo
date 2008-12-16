function R = rndalpha(varargin)
%RNDALPHA Random matrices from a symmetric alpha-stable distribution
% 
% CALL:  R = rndalpha(a,sz);
%
%        R = matrix of random numbers, 
%        a = the parameter  alpha, 
%       sz = size(R)    (Default size(a))
%            sz is a comma separated list or a vector 
%            giving the size of R (see zeros for options).
%
% The characteristic function of a symmetric alpha-stable distribution is
%   h(t) = exp(-abs(t)^a)
% a = 1 gives the Cauchy distribution and a = 2 the Gaussian distribution.
%
% Example:
%   R=rndalpha(0.5,1,100);
%   plot(R,'.')
%   R=rndalpha(2,1,100);
%   plotnorm(R)
% 
% See also  zeros

% Reference:
% Samordnitsky & Taqqu (1994) "Non-Gaussian and Stable Processes"
% Chapman & Hall

% Tested on: Matlab 5.3
% History: 
% Revised by jr 22.12.1999
% revised ms 15.06.2000
% - updated header info
% - changed name to rndalpha (from alpharnd)
% - rewrote code to work with matrices
% revised pab 23.10.2000
%  - added comnsize, nargchk
%  - added greater flexibility on the sizing of R

error(nargchk(1,inf,nargin))
Np = 1;
options = []; %struct; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
% if numel(options)>1
%   error('Multidimensional struct of distribution parameter not allowed!')
% end

a = params{1};
if isempty(rndsize)
  csize=size(a);
else
  [csize] = comnsize(a,zeros(rndsize{:}));
  if any(isnan(csize))
    error('a  must be  scalar or conform to the size information given.');
  end
end


Y1=-log(rand(csize));
Y2=pi*(rand(csize)-0.5);

R=sin(a.*Y2)./cos(Y2).^(1./a).*(cos((1-a).*Y2)./Y1).^((1-a)./a);
