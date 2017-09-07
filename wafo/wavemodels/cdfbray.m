function [F,Flo,Fhi] = cdfbray(x,varargin)
%CDFBRAY Beta Rayleigh CDF of wave heights 
%       /h              
%   F = | 2*(a+b-1)!/((a-1)! * (b-1)!)*x^(2a-1)*(1-(x/c)^2)^(b-1)/c^(2a) dx
%       /0                
%
%  CALL:  F = cdfbray(h,a,b,c) 
%
%       F = cdf
%       h = waveheigth (0 <= h <= c)
%       a = abs(k1*(k2-k1)/(k1^2-k2)) 
%       b = abs(1-k1)*(k2-k1)/(k1^2-k2)) 
%       c = Hb, breaking wave height approximated by water depth, d.
% where
%      k1 = E(H^2)/Hb^2
%      k2 = E(H^4)/Hb^4
%  E(H^2) = .5*exp(0.00272*(d/g*Tp^2)^(-0.834))*Hm0^2
%  E(H^4) = .5*exp(0.00046*(d/g*Tp^2)^(-1.208))*Hm0^2
%     Hm0 = significant waveheight
%     Tp  = modal period of wave spectrum
%
%    The size of F is the common size of H, A, B and C.  A scalar input   
%    functions as a constant matrix of the same size as the other input.
%
% Example: % Compare with rayleigh distribution
%  Hm0 = 7;Tp = 11;d = 50; g = gravity;
%  k1  = .5*exp(0.00272*(d/g*Tp^2)^(-0.834))*Hm0^2/d^2;
%  k2  = .5*exp(0.00046*(d/g*Tp^2)^(-1.208))*Hm0^2/d^4;
%  a   = abs(k1*(k2-k1)/(k1^2-k2)); 
%  b   = abs((1-k1)*(k2-k1)/(k1^2-k2));
%  h   = linspace(0,2*Hm0)';
% semilogy(h,1-cdfbray(h,a,b,d),'r',h,1-cdfray(h,Hm0/2))
%
% See also  cdfbeta

% 
%   Reference:
%       Michel K. Ochi (1998),
%      "OCEAN WAVES, The stochastic approach",
%       OCEAN TECHNOLOGY series 6, Cambridge, pp 279. (pd of peaks to trough) 

% tested on: matlab 5.x
% History:
% Revised pab 31.03.2001 
%  added example
% revised pab 14.10.1999
% updated help header
%  Per A. Brodtkorb 21.02.99
%error(nargchk(4,12,nargin))
narginchk(4,12)
options = struct('alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});
a(a<=0) = nan;
b(b<=0) = nan;
c(c<=0) = nan;
x(x<0)  = 0;
try
  xn = (x./c).^2;
  xn(xn>1) = 1;
  if nargout>1
    phat.params(end) = [];
    phat.covariance(:,3) = [];
    phat.covariance(3,:) = [];
    [F,Flo,Fhi] = cdfbeta(xn,phat,options);
  else
    F = cdfbeta(xn,a,b,options);
  end
catch
    error('h, a, b and c must be of common size or scalar.');
end






