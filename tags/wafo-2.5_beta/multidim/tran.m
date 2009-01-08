function [Hw, Gwt, kw,Hwt,ee]=tran(w,theta,pos,def,h,g,rho,bet,igam,thx,thy,kw)
%TRAN Computes transfer functions based on linear wave theory
%     of the system with input surface elevation, 
%     eta(x0,y0,t) = exp(i*(kx*x0+ky*y0-w*t)), 
%     and output Y determined by type and pos. 
%
% CALL:  [Hw Gwt kw Hwt] = tran(w,theta,pos,type,h,g,rho,bet,igam,thx,thy,kw);
%
%   Hwt = Hw(ones(Nt,1),:).*Gwt matrix of transfer functions values as function of 
%         w (columns) and theta (rows)                   size Nt x Nf 
%   Hw  = a function of frequency only (not direction)   size  1 x Nf
%   Gwt = a function of frequency and direction          size Nt x Nf
%     w = vector of angular frequencies in Rad/sec. Length Nf
% theta = vector of directions in radians           Length Nt   (default 0)
%         ( theta = 0 -> positive x axis theta = pi/2 -> positive y axis)
%   pos = [x,y,z] = vector giving coordinate position relative to [x0 y0 z0] (default [0,0,0])
%  type = Number or string defining the sensortype or transfer function in output. 
%         1,  'n'    : Surface elevation              (n=Eta)     (default)
%         2,  'n_t'  : Vertical surface velocity
%         3,  'n_tt' : Vertical surface acceleration
%         4,  'n_x'  : Surface slope in x-direction 
%         5,  'n_y'  : Surface slope in y-direction
%         6,  'n_xx' : Surface curvature in x-direction
%         7,  'n_yy' : Surface curvature in y-direction
%         8,  'n_xy' : Surface curvature in xy-direction
%         9,  'P'    : Pressure fluctuation about static MWL pressure 
%         10, 'U'    : Water particle velocity in x-direction
%         11, 'V'    : Water particle velocity in y-direction
%         12, 'W'    : Water particle velocity in z-direction
%         13, 'U_t'  : Water particle acceleration in x-direction
%         14, 'V_t'  : Water particle acceleration in y-direction
%         15, 'W_t'  : Water particle acceleration in z-direction
%         16, 'X_p'  : Water particle displacement in x-direction from its mean position
%         17, 'Y_p'  : Water particle displacement in y-direction from its mean position
%         18, 'Z_p'  : Water particle displacement in z-direction from its mean position
%         19, 'Dummy': Transfer function is zero
%    h  = water depth      (default inf) 
%    g  = acceleration of gravity (default see gravity)
%   rho = water density    (default see wdensity)
%   bet = 1, theta given in terms of directions toward which waves travel (default)
%        -1, theta given in terms of directions from which waves come 
%  igam = 1, if z is measured positive upward from mean water level (default)
%         2, if z is measured positive downward from mean water level
%         3, if z is measured positive upward from sea floor
%   thx = angle clockwise from true north to positive x-axis in degrees
%         (default 90)
%   thy = angle clockwise from true north to positive y-axis in degrees
%         (default 0)
%    kw = vector of wave numbers corresponding to angular frequencies, w
%         (default calculated with w2k)
%
% Example:
%   N=5000;dt=0.4;f0=0.1;th0=0;h=50; w0 = 2*pi*f0;
%  t = linspace(0,1500,N)';
%  types = [1 4 5];
%  bfs   = [1 0 0];
%  pos   = zeros(3,3);
%  eta0 = exp(-i*w0*t);
%  [Hw Gwt] = tran(w0,th0,pos(1,:),'n',h); 
%  eta = real(Hw*Gwt*eta0)+0.001*rand(N,1);
%  [Hw Gwt] = tran(w0,th0,pos(2,:),'n_x',h); 
%  n_x = real(Hw*Gwt*eta0)+0.001*rand(N,1);
%  [Hw Gwt] = tran(w0,th0,pos(3,:),'n_y',h); 
%  n_y = real(Hw*Gwt*eta0)+0.001*rand(N,1);
%
%  S = dat2dspec([t eta n_x n_y],[pos types',bfs'],h,256,101);
%  plotspec(S)
%
% See also  dat2dspec, sensortype, sensortypeid, wdensity, gravity

% Reference
% Young I.R. (1994)
% "On the measurement of directional spectra",
% Applied Ocean Research, Vol 16, pp 283-294

% History
% revised pab jan2005 
% -moved the code changing sensortypename to sensortypeid into a separate
% function sensortypeid.m
% revised pab 11.10.2002
% - fixed some bugs in transfer functions
% - reordered the transfer functions
% - 
% revised pab 07.11.2001
% -Made a fix for w=0 where the transfer functions sometimes is singular.
%  e.g., Hw  = w.*ratio(zk,hk,-1,-1); return NaN for zk=0,hk=0 and w=0.
%  This fix is strictly not correct for all cases, but suffices in most cases.
% -moved tran as a local function to a function in \wafo\multidim\
% revised pab 14.06.2000
% - updated documentation
% - clearified directions
% - added def for surface curvature
% revised pab 06.01.2000
%  - cleaned up documentation
%  - added default values
%  - speeded up computations
% based on transp by L. Borgman

% Default values
%~~~~~~~~~~~~~~~
if nargin<2||isempty(theta), theta = 0; end
if nargin<3||isempty(pos),   pos   = [0 0 0]; end
if nargin<4||isempty(def),   def   = 1; end  
if nargin<5||isempty(h),     h     = inf; end
if nargin<6||isempty(g),     g     = gravity; end % accelleration of gravity
if nargin<7||isempty(rho),   rho   = wdensity; end % water density given in kg/m^3
if nargin<8||isempty(bet),   bet   = 1;  end
if nargin<9||isempty(igam),  igam  = 1;  end
if nargin<10||isempty(thx),  thx   = 90; end
if nargin<11||isempty(thy),  thy   = 0;  end
if nargin<12||isempty(kw),   kw    = w2k(w,0,h); end % find the wave number as function of angular frequency

if ischar(def)
  def1 = sensortypeid(def);
  if isempty(def1)||length(def1)>1
    error('Invalid input def = %s',def)
  end
  def = def1;
end

% make sure they have the correct orientation
theta = theta(:);
kw    = kw(:).';
w     = w(:).';


% Compute various fixed arrays
% ----------------------------
Nt    = length(theta);
Nf    = length(w);
Oneth = ones(Nt,1);
Onef  = ones(1,Nf);



% convert to from angle in degrees to radians
thxr = thx*pi/180;
thyr = thy*pi/180;

cthx = bet*cos(theta-thxr+pi/2);
%cthy = cos(theta-thyr-pi/2);
cthy = bet*sin(theta-thyr);

% Compute location complex exponential
% ------------------------------------

x=pos(1);y=pos(2);z=pos(3);
%arg = bet*(x*cthx+y*cthy)*kw;
%arg = arg(:,Onef).*kw(Oneth,:);

% old call
%ee=cos(arg)-(sqrt(-1))*sin(arg); % exp(-i*k(w)*(x*cos(theta)+y*sin(theta)) size Nt X Nf
%new call
ee= exp((sqrt(-1)*(x*cthx+y*cthy))*kw);   % exp(i*k(w)*(x*cos(theta)+y*sin(theta)) size Nt X Nf


if def > 8
   % Set up z and h arguments
   % ------------------------ 
   hk = kw*h;
   switch igam
     case 1,   zk = kw*(h+z); % z measured positive upward from mean water level (default)
     case 2,   zk = kw*(h-z); % z measured positive downward from mean water level
     otherwise,zk = kw*z;     % z measured positive upward from sea floor
   end
end

% Start computation of complex transformations
% ---------------------------------------------
switch def
  case 1,                                    % n   = Eta = wave profile
    Hw  = Onef;
    Gwt = ee;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %             Vertical surface velocity and acceleration
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 2,                                    % n_t = Eta_t 
    Hw  = w;
    Gwt = -sqrt(-1)*ee;
  case 3,                                     % n_tt = Eta_tt
    Hw  = w.^2;
    Gwt = -ee;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %             Surface slopes
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 4,                                    % n_x = Eta_x = x-slope
    Hw  = kw;
    Gwt = sqrt(-1)*cthx(:,Onef).*ee;
  case 5,                                    % n_y = Eta_y = y-slope
    Hw  = kw;
    Gwt = sqrt(-1)*cthy(:,Onef).*ee;
    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %             Surface curvatures
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 6,                                   % n_xx = Eta_xx = Surface curvature (x-dir)
    Hw = kw.^2;
    Gwt = -cthx(:,Onef).^2.*ee; 
  case 7,                                   % n_yy = Eta_yy = Surface curvature (y-dir)
    Hw = kw.^2;
    Gwt = -cthy(:,Onef).^2.*ee; 
  case 8,                                   % n_xy = Eta_xy = Surface curvature (xy-dir)
    Hw = kw.^2;
    Gwt = -cthx(:,Onef).*cthy(:,Onef).*ee; 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %            Pressure
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
  case 9, % pressure fluctuations
    Hw  = rho*g*ratio(zk,hk,1,1);       % cosh(zk)/cosh(hk)
    Gwt = ee;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %             water particle velocities
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 10,                                    % U = x-velocity
    Hw  = w.*ratio(zk,hk,1,-1); % w*cosh(zk)/sinh(hk)
    Gwt = cthx(:,Onef).*ee; % cos(theta)*ee 
  case 11,                                    % V = y-velocity
    Hw  = (w.*ratio(zk,hk,1,-1)); % w*cosh(zk)/sinh(hk)
    Gwt = cthy(:,Onef).*ee; % sin(theta)*ee 
  case 12,                                    % W = z-velocity
    Hw  = w.*ratio(zk,hk,-1,-1); % w*sinh(zk)/sinh(hk)
    Gwt = -sqrt(-1)*ee;                        %-?         
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %             water particle acceleration
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  case 13 ,                                   % U_t = x-acceleration
    Hw  = (w.^2).*ratio(zk,hk,1,-1);% w^2*cosh(zk)/sinh(hk)
    Gwt = -sqrt(-1)*cthx(:,Onef).*ee;                         %-?
  case 14  ,                                  % V_t = y-acceleration
    Hw  = (w.^2).*ratio(zk,hk,1,-1);% w^2*cosh(zk)/sinh(hk)
    Gwt = -sqrt(-1)*cthy(:,Onef).*ee;                         %-?
  case 15,                                    % W_t = z-acceleration
    Hw  = (w.^2).*ratio(zk,hk,-1,-1);% w*sinh(zk)/sinh(hk)
    Gwt = -ee;
    
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %             water particle displacement
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  case 16,                                    % X_p = x-displacement
    Hw  = ratio(zk,hk,1,-1);            % cosh(zk)./sinh(hk)
    Gwt = sqrt(-1)*cthx(:,Onef).*ee;
  case 17,                                    % Y_p = y-displacement
    Hw  = ratio(zk,hk,1,-1);            % cosh(zk)./sinh(hk) 
    Gwt = sqrt(-1)*cthy(:,Onef).*ee; 
  case 18,                                    % Z_p = z-displacement
    Hw  = ratio(zk,hk,-1,-1);            % sinh(zk)./sinh(hk) 
    Gwt = ee;
  otherwise
    error('Unknown option 1<=def<=18')
 end
 
 % New call to avoid singularities. pab 07.11.200
 % Set Hw to 0 for expressions w*ratio(z*k,h*k,1,-1)= 0*inf
Hw(isnan(Hw) | isinf(Hw)) = 0;

sgn = sign(Hw);
k0 = find(sgn<0);
if any(k0) % make sure Hw>=0 ie. transfer negative signs to Gwt
  Gwt(:,k0) = -Gwt(:,k0);
  Hw(k0)    = -Hw(k0);
  
end
if igam==2, 
  %pab 09 Oct.2002: bug fix
  %  Changing igam by 2 should affect the directional result in the same way that changing eta by -eta!
  Gwt = -Gwt;
end


if nargout>3
 Hwt=Hw(Oneth,:).*Gwt;
end
%if Hw is wanted in size Nt x Nf uncomment this!!!
%   Hw=Hw(Oneth,:); 

return

