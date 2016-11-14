function D = testsurf(N,dt,amp,f0,thet0,x,y,d,g,thetx,thety)
%TESTSURF creates a test case for a surface elevation measurement
% 
% CALL  D = testsurf(N,dt,amp,f0,thet0,x,y,h,g,thetx,thety);
%
%   D     = matrix containing column vectors of time series
%   N     = number of time steps
%   dt    = time increment
%   amp   = amplitude
%   f0    = primary frequency in Hz
%   thet0 = primary direction in degrees toward which the waves travels
%           (0 = North, 90 = East,...etc)
%   x     = vector of x-coordinates length M
%   y     = vector of y-coordinates length M
%   h     = water depth               (default infinity)
%   g     = acceleration of gravity   (default see gravity)
%   thetx = angle in degrees clockwise from north to the + x-axis
%                   (default 90)
%   thety = angle in degrees clockwise from north to the + y-axis
%                   (default 0)
%
% CREATE A TEST CASE:
%  eta = amp*cos(k*x*cos(th0)+k*y*sin(th0)-2*pi*f0*t);
%     with cos(th0) = cos(thet0-thetx);
%          sin(th0) = cos(thet0-thety);
% Example:
%   N=5000;dt=0.4;f0=0.1;th0=0;h=50;xypos = [0 0 0 1 1;0 40 0 1 1; 20 20 0 1 1];
%   D = testsurf(N,dt,3,f0,th0,xypos(:,1),xypos(:,2),h);
%   S = dat2dspec(D,xypos,h);
%   plotspec(S);
%
%   close all;
%
% See also  testbuoy, dat2dspec

% History:
% revised pab 14.06.2000
% - updated documentation
% revised pab 06.01.2000
%  - updated documentation
%  - added default values
% based on testbuoy by  L. Borgman

% NOTE : Directions towards which wave travels are correct! PAB Feb 2007

if nargin<8||isempty(d)
 d=inf;
end
if nargin<9||isempty(g)
 g = gravity;
end
if nargin<10||isempty(thetx)
 thetx = 90;
end
if nargin<11||isempty(thety)
 thety = 0;
end

thet0r = thet0*pi/180;
thetxr = thetx*pi/180;
thetyr = thety*pi/180;
%[0 40 20],[0 0 40]

X = ones(N,1)*x(:)';
Y = ones(N,1)*y(:)';
t = (0:N-1)'*dt;
T = t*ones(1,length(x(:)));

% Compute wave number
% -------------------
k = w2k(2*pi*f0,0,d,g);

% Compute test time series
% ------------------------
%cos(thet0r-thetxr)
%cos(thet0r-thetyr)

if 1, % orientation clockwise th0 = 0 = north th0=90= east
  D = amp*cos(k*X*cos(thet0r-thetxr)+k*Y*cos(thet0r-thetyr)-2*pi*f0*T);
else  % orientation anti clockwise  th0=0 = east th0=90 = north
  D = amp*cos(k*X*cos(thet0r-thetxr+pi/2)+k*Y*sin(thet0r-thetyr)-2*pi*f0*T);
end
% Add in some noise
% -----------------
D = D+.01*sqrt(((amp^2)/2)/1000)*randn(size(D));
D = [t D];
