function D = testbuoy(N,dt,amp,f0,thet0,x,y,d,g,thetx,thety)
%TESTBUOY Creates a test case for a buoy measurement
% 
% CALL  D = testbuoy(N,dt,amp,f0,thet0,x,y,h,g,thetx,thety);
%
%   D     = matrix of time series 
%           for M=1:
%           c1: time c2: eta c3: x-slope c4: y-slope
%   N     = number of time steps
%   dt    = time increment
%   amp   = amplitude
%   f0    = primary frequency in Hz
%   thet0 = primary direction in degrees toward which the waves travels
%           (0 = North, 90 = East)
%   x     = vector of x-coordinates length M
%   y     = vector of y-coordinates length M
%   h     = water depth               (default infinity)
%   g     = accelleration of gravity  (default see gravity)
%   thetx = angle in degrees clockwise from north to the + x-axis
%                   (default 90)
%   thety = angle in degrees clockwise from north to the + y-axis
%                   (default 0)
%
% CREATE BUOY TEST CASE:
%  eta = amp*cos(k*x*cos(th0)+k*y*sin(th0)-2*pi*f0*t);
%  x-slope = -amp*k*cos(th0)*sin(k*x*cos(th0)+k*y*sin(th0)-2*pi*f0*t);
%  y-slope = -amp*k*sin(th0)*sin(k*x*cos(th0)+k*y*sin(th0)-2*pi*f0*t);
%     with cos(th0) = cos(thet0-thetx);
%          sin(th0) = cos(thet0-thety);
%
% Example:
%   N=5000;dt=0.4;f0=0.1;th0=0;h=50;xypos=[0 0 0 1 1;0 0 0 4 0; 0 0 0 5 0];
%   D = testbuoy(N,dt,3,f0,th0,0,0,h);
%   S = dat2dspec(D,xypos,h);
%   plotspec(S);
%
%   close all;
%
% See also  testsurf, dat2dspec

% History:
% revised pab 14.06.2000
% updated documentation
% revised pab 06.01.2000
%  - updated documentation
%  - added default values
%  - corrected the D matrix + added sampling times
% by  L. Borgman

if nargin<8||isempty(d),      d     = inf;    end
if nargin<9||isempty(g),      g     = gravity;end
if nargin<10||isempty(thetx), thetx = 90;     end
if nargin<11||isempty(thety), thety = 0;      end

thet0r=thet0*pi/180;
thetxr=thetx*pi/180;
thetyr=thety*pi/180;

X=ones(N,1)*x(:)';
Y=ones(N,1)*y(:)';
t=(0:N-1).'*dt;
T=t*ones(1,size(x,1));

% Compute wave number
% -------------------
k=w2k(2*pi*f0,0,d,g);

% Compute test time series
% ------------------------
mx=k*cos(thet0r-thetxr);
my=k*cos(thet0r-thetyr);
D=amp*cos(k*X*cos(thet0r-thetxr)+k*Y*cos(thet0r-thetyr)-2*pi*f0*T);
D=[D, -mx*amp*sin(k*X*cos(thet0r-thetxr)+k*Y*cos(thet0r-thetyr)-2*pi*f0*T)];
D=[D, -my*amp*sin(k*X*cos(thet0r-thetxr)+k*Y*cos(thet0r-thetyr)-2*pi*f0*T)];

% Add in some noise
% -----------------
D=D+.01*sqrt(((amp^2)/2)/100)*randn(size(D));
D=[T,D];
