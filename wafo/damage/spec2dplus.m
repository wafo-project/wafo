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
%  dplus = spec2dplus(S,bet)

% Tested on: Matlab 5.3
% History: 
% By jr 12.01.2000
  
lam = spec2mom(S); l0 = sqrt(lam(1));, l2 = sqrt(lam(2));

d = (2.^(bet/2)).*(l0.^(bet-1)).*gamma(bet/2+1)*l2/2/pi;
