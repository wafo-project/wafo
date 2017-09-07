function rho = wdensity(S,T,P)
%WDENSITY  Returns the water density 
%
%  CALL:  rho = wdensity(S,T,P);
%
%   rho =  water density [kg/m^3]
%   S   = salinity       [psu]       (Default 35)
%   T   = temperature    [degrees C] (default 4)
%   P   = pressure       [db]        (default 0) 
%
% WDENSITY estimates the water density as function of salinity, temperature
% and pressure. The formulae is valid in the region:
%
%  0 < S [psu] < 42, -2 < T [C] < 40 and 0 < P [db] < 10000
%
% Example:
%  S = linspace(20,40)'; T = linspace(4,20)';
%  [S1 T1]=meshgrid(S,T);
%  sc = contour(S,T,wdensity(S1,T1)); clabel(sc)
%  xlabel('Salinity'),ylabel('Temperature')

% REFERENCES:
%  Fofonoff, P. and Millard, R.C. Jr (1983)
%  Algorithms for computation of fundamental properties of 
%  seawater, Unesco Tech. Pap. in Mar. Sci., No. 44, 53, pp 17-18
%
%  Millero, F.J & Poisson, A. (1981)
%  International one-atmosphere equation of state for seawater.
%  Deep-Sea Research, Vol. 28A, No.6, pp 625-629. 
%
%  BIGG P.H.,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537.
%

% Copyright (C) 2000  Per A. Brodtkorb
% 
%  This file, WDENSITY.M, is part of WAFO.
% 
%     WDENSITY is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     WDENSITY is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%error(nargchk(0,3,nargin))
narginchk(0,3)
if nargin<1||isempty(S), 
  S = 35; 
end
if nargin<2||isempty(T),
  T = 4; 
end
if nargin<3||isempty(P),
  P = 0; 
end
if any(S <0 | 42 < S)
  warning('WAFO:WDENSITY','Salinity is outside valid boundaries [0 42] \n The validity of the formulae is questionable') 
end
if any(T <-2 | 40 < T)
  warning('WAFO:WDENSITY','Temperature is outside valid boundaries [-2 40] \n The validity of the formulae is questionable') 
end
if any(P <0 | 10000 < P)
  warning('WAFO:WDENSITY','Pressure is outside valid boundaries [0 10000] \n The validity of the formulae is questionable') 
end


[iscsize, S,T,P] = iscomnsize(S,T,P);
if ~iscsize>0,
  error('S, T and P must be scalar or of common size')
end


% Density of fresh water at atmospheric pressure
%   BIGG P.H.,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537.
%-------------------------------------------------
a  = [6.536332e-9, -1.120083e-6, 1.001685e-4, -9.095290e-3, 6.793952e-2 999.842594];
rho0 = polyval(a,T); % unit: [kg/m^3]

% Density of salt water at atmospheric pressure:
%------------------------------------------------
acof = [ 5.3875e-9, -8.2467e-7,  7.6438e-5, -4.0899e-3, 8.24493e-1];
bcof = [ -1.6546e-6, 1.0227e-4, -5.72466e-3];
  

A = polyval(acof,T);
B = polyval(bcof,T);
C = 4.8314e-4;
rho = rho0 + (A + B.*sqrt(S) + C*S).*S;


k2 = find(P~=0);
if any(k2),
  Pk  = P(k2);
  K   = seck(S(k2),T(k2),Pk);
  Pk  = Pk/10;  % convert from db to atmospheric pressure units
  rho(k2) = rho(k2)./(1-Pk./K);
end
return


function K = seck(S,T,P)
%SECK  Secant bulk modulus (K) of sea water
%
% CALL:  K = seck(S,T,P);
%
%   S = salinity    [psu]
%   T = temperature [degrees C]
%   P = pressure    [db]
%
%   K = Secant Bulk Modulus  [bars]
% 
%    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980. 
%    UNESCO polynomial implementation.

% REFERENCES:
%  Fofonoff, P. and Millard, R.C. Jr (1983)
%  Algorithms for computation of fundamental properties of 
%  seawater, Unesco Tech. Pap. in Mar. Sci., No. 44, 53, pp 17-18

%error(nargchk(3,3,nargin))
narginchk(3,3)

% Fresh water terms of K at atmosphere pressure.
%------------------------------------------------
e = [-5.155288E-5, 1.360477E-2, -2.327105, 148.4206,19652.21]; 
K  = polyval(e,T);   % Eq. 19

% Salt water terms at atmosphere pressure
%---------------------------------------------
SR = sqrt(S);
f = [-6.1670E-5, 1.09987E-2,-0.603459,54.6746];
g = [-5.3009E-4, 1.6483E-2,7.944E-2];
K = K + ( polyval(f,T) + polyval(g,T).*SR ).*S;    % eq. 16

%
%---------------------------------------------
k2 = find(P~=0);
if any(k2),
  Tk = T(k2);
  h  = [-5.77905E-7, 1.16092E-4, 1.43713E-3, 3.239908];
  j0 = 1.91075E-4;
  icoef =[ -1.6078E-6, -1.0981E-5  2.2838E-3];
  A  = polyval(h,Tk) + (polyval(icoef,Tk)+ j0*SR(k2)).*S(k2); 
  
  kcoef = [ 5.2787E-8, -6.12293E-6,8.50935E-5];  
  m = [ 9.1697E-10, 2.0816E-8, -9.9348E-7];
  B = polyval(kcoef,Tk) + polyval(m,Tk).*S(k2);  
  
  Pk = P(k2)/10;  %convert from decibars to atmospheric pressure units (bars)
  K(k2) = K(k2) + (A + B.*Pk).*Pk;  % eq. 15
end
return









