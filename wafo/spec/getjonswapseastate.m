function seastate = getjonswapseastate(u10,fetch,method)
% GETJONSWAPSEASTATE Estimate seastate from windspeed and fetch
%
%  CALL: seastate = getjonswapseastate(U10,fetch,method)
%
% seastate = [Hm0,Tp, gamma sa sb A] where
%            Hm0   : significant wave height [m]
%            Tp    : peak period [s]
%            gamma : jonswap peak enhancement factor.
%            sa,sb : jonswap spectral width parameters. 
%            A     : jonswap alpha, normalization factor.
% U10    = windspeed at 10 m above mean water surface [m/s]
% fetch  = fetch [m] (default 150000)
% method = 'hasselman73' seastate according to Hasselman et. al. 1973
%          'hasselman76' seastate according to Hasselman et. al. 1976
%          'lewis'       seastate according to Lewis and Allos 1990 (default)
% 
%Example
% fetch = 10000;
% u10   = 10;
% ss = getjonswapseastate(u10,fetch)
% S  = jonswap([],ss);
% plotspec(S)
% spec2char(S,{'hm0','Tp'})
%
% See also jonswap


% References
% Lewis, A. W. and Allos, R.N. (1990)
% JONSWAP's parameters: sorting out the inconscistencies.
% Ocean Engng, Vol 17, No 4, pp 409-415
%
% Hasselmann et al. (1973)
% Measurements of Wind-Wave Growth and Swell Decay during the Joint
% North Sea Project (JONSWAP). 
% Ergansungsheft, Reihe A(8), Nr. 12, Deutschen Hydrografischen Zeitschrift.
%
% Hasselmann et al. (1976)
% A parametric wave prediction model.
% J. phys. oceanogr. Vol 6, pp 200-228

%History
% by pab 6 jan2006
%

error(nargchk(1,3,nargin))
if nargin < 2 || isempty(fetch)
  fetch = 150000;
end
if nargin<3 || isempty(method)
  method = 'lewis';
end
g    = gravity;

% The following formulas are from Lewis and Allos 1990:

zeta    = g*fetch./(u10^2); % dimensionless fetch, Table 1
%zeta = min(zeta, 2.414655013429281e+004);  

if strncmpi(method,'hasselman',1)
  if method(end)=='3',
    % Hasselman et.al (1973)
    A       = 0.076*zeta^(-0.22);
    ny      = 3.5*zeta^(-0.33);               % dimensionless peakfrequency, Table 1
    epsilon1 = 9.91e-8*zeta.^1.1;             % dimensionless surface variance, Table 1
  else
    % Hasselman et.al (1976)
    A       = 0.0662*zeta^(-0.2);
    ny      = 2.84*zeta^(-0.3);               % dimensionless peakfrequency, Table 1
    epsilon1 = 1.6e-7*zeta;                   % dimensionless surface variance, Eq.4
  end
  sa      = 0.07;
  sb      = 0.09;
  gam     = 3.3;
else
  A       = 0.074*zeta^(-0.22);               % Eq. 10
  ny      = 3.57*zeta^(-0.33);                % dimensionless peakfrequency, Eq. 11
  epsilon1 = 3.512e-4*A*ny^(-4)*zeta^(-0.1);  % dimensionless surface variance, Eq.12
  sa      = 0.05468*ny^(-0.32);               % Eq. 13
  sb      = 0.078314*ny^(-0.16);              % Eq. 14
  gam     = max(17.54*zeta^(-0.28384),1);     % Eq. 15
end
Tp      = u10/(ny*g);                          % Table 1
Hm0     = 4*sqrt(epsilon1)*u10^2/g;            % Table 1
seastate = [Hm0, Tp, gam, sa, sb, A];