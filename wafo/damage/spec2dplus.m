function d = spec2dplus(S,bet)
% SPEC2DPLUS Calculates an upper bound of the damage intensity explicitly
%
% CALL:  d = spec2dplus(S,bet)
%
%        d = an upper bound of the damage intensity
%      
%        S = a spectral density (structure array)
%      bet = the parameter beta
%
% The upper bound of the damage intensity can be calculated explicitly and is
% given by 
% 
%       d^+ = 2^{3*beta/2}*sqrt(lam2)*(lam0)^{(beta-1)/2}*gamma(1+beta/2)/(2*pi)
%
% Example
%  S = oscspec; bet = 3:0.2:5;
%  dplus = spec2dplus(S,bet);
%  assert(dplus, [ 2.32262632058483, 2.93360907105837,3.72269540324765,...
%                  4.74521181972958, 6.07455506745236, 7.80833828497044,...
%                 10.07670105706864, 13.05356835488749, 16.97193573771180,...
%                 22.14466211791423, 28.99281030590514], 1e-10);

% Tested on: Matlab 5.3
% History: 
% By jr 12.01.2000
  
lam = spec2mom(S); l0 = sqrt(lam(1));, l2 = sqrt(lam(2));

d = (2.^(bet/2)).*(l0.^(bet-1)).*gamma(bet/2+1)*l2/2/pi;
