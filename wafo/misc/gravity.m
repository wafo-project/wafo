function g=gravity(phi)
% GRAVITY  returns the constant acceleration of gravity 
%
%  CALL:  g = gravity(phi);
%
%   g   =  acceleration of gravity in m/s^2 
%   phi =  latitude in degrees (default 45 degrees) 
%
%  GRAVITY calculates the acceleration of gravity 
%  using the international gravitational formulae:
%
%   g = 9.78049*(1+0.0052884*sin(phi).^2-0.0000059*sin(2*phi).^2);
%
%  Edit GRAVITY.M to change default value for PHI.
%
% Example:%  
%    g = gravity
% See also wdensity

% References:
% Irgens, Fridtjov (1987)
% "Formelsamling i mekanikk: 
% statikk, fasthetslære, dynamikk fluidmekanikk"
% tapir forlag, University of Trondheim , ISBN 82-519-0786-1, pp 19

% history
% revised pab 06.01.2000
%  - added phi and international gravitational formulae.
%  - updated documentation

if nargin<1 ||isempty(phi)
  phi=45;
end
phi=phi*pi/180; % change from degrees to radians

g=9.78049*(1+0.0052884*sin(phi).^2-0.0000059*sin(2*phi).^2); 
%g=9.80665000000;
