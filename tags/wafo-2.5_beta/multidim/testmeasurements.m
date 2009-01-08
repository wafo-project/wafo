function D = testmeasurements(pos,type,thet0,f0,N,dt,amp,d,g,thetx,thety)
%TESTMEASUREMENTS Creates a test case for measurement time series
% 
% CALL:  D = testmeasurements(pos,type,thet0,f0,N,dt,amp,h,g,thetx,thety);
%
%   D     = matrix containing column vectors of time series, size N x M+1
%   pos   = coordinate position of the sensors in each row, size M x 3.
%   type  = vector of sensortypes given as integers see tran for options, length M
%   thet0 = primary direction in degrees toward which the waves travels
%           (0 = East, 90 = North,...etc)  (default 0)
%   f0    = primary frequency in Hz        (default 0.1 Hz)
%   N     = number of time steps           (default 5000)   
%   dt    = time increment                 (default 0.5)
%   amp   = amplitude                      (default 1)
%   h     = water depth                    (default infinity)
%   g     = acceleration of gravity        (default see gravity)
%   thetx = angle in degrees clockwise from north to the + x-axis
%                                          (default 90)
%   thety = angle in degrees clockwise from north to the + y-axis
%                                          (default 0)
%
% CREATE A TEST CASE:
%  eta = amp*cos(k*x*cos(th0)+k*y*sin(th0)-2*pi*f0*t);
%     with cos(th0) = cos(thet0-thetx);
%          sin(th0) = cos(thet0-thety);
% Example:
%   type = [1 1 1]; bfs = ones(1,3);h=inf;
%   th0  = 90;
%   pos = [0 0 0;0 40 0; 20 20 0];
%   D = testmeasurements(pos,type,th0);
%   S = dat2dspec(D,[pos type' bfs'],h);
%   plotspec(S)
%
% See also  testbuoy, dat2dspec, tran

% History:
% by pab 14.10.2002


error(nargchk(2,11,nargin))

[M ] = size(pos,1);
if M~=length(type)
    error('size(pos,1) must be equal to length(type)')
end

if nargin<3||isempty(thet0),  thet0 = 0;end
if nargin<4||isempty(f0),     f0    = 0.1;end
if nargin<5||isempty(N),      N     = 5000;end
if nargin<6||isempty(dt),     dt    = 0.5;end
if nargin<7||isempty(amp),    amp   = 1;end
if nargin<8||isempty(d),      d     = inf;end
if nargin<9||isempty(g),      g     = gravity;end
if nargin<10||isempty(thetx), thetx = 90;end
if nargin<11||isempty(thety), thety = 0;end


%convert from degrees to radians
thet0r=thet0*pi/180;




% Compute wave number
% -------------------
w0 = 2*pi*f0;
kw = w2k(w0,0,d,g);

t=(0:N-1)'*dt;
eta0 = amp*exp(-i*w0*t);

D = zeros(N,M);

for ix=1:M
     [Hw, Gwt] = tran(w0,thet0r,pos(ix,:),type(ix),d,g,[],[],[],thetx,thety,kw); 
     D(:,ix) = real(Hw*Gwt*eta0)+0.0001*sqrt(amp^2)*rand(N,1);
end 

D = [t D];


