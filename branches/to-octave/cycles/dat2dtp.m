function [dtp,u,tp] = dat2dtp(param,x,h,ddef)
%DAT2DTP The sequence of discretized turning points from a signal.
%
% CALL:  [dtp,u,tp] = dat2dtp(param,x,h,ddef);
%               dtp = dat2dtp(param,x);
%     
% Input:
% param = Parameter vector, [a b n], defines the discretization.
% x     = Data, one column with values OR                   [N,1]
%               two columns with times and values.          [N,2]
% h     = Threshold for hysteresis filter (see rfcfilter). 
%           (Optional, Default='smallest discretization step')
% ddef  = 1 causes peaks to be projected upwards and troughs 
%           downwards to the closest discrete level (default).
%       = 0 causes peaks and troughs to be projected
%           the closest discrete level.
%       =-1 causes peaks to be projected downwards and the 
%           troughs upwards to the closest discrete level.
%
% Output:
% dtp   = Discretized turning points taking values 1,2,...,n.  
%                                                   [N,1] or [N,2]
% u     = Discretization levels.                             [1,n]
% tp    = Discretized turning points taking values u(1),u(2),...,u(n).
%                                                   [N,1] or [N,2]
%
% Calculates the sequence of discretized turning points and optionally 
% removes small oscillations from data x by rainflow filtering.
%
% Example:
%   x = load('sea.dat'); x=x(1:200,:);     % Load data
%   [dtp,u,tp] = dat2dtp([-2 2 32],x,0.5); % Discrete TP & rainflow filter 0.5
%   plot(x(:,1),x(:,2),tp(:,1),tp(:,2))
%
% See also  dat2tp, tp2cc, dtp2rfm, rfcfilter.

% Copyright (c) 1999 by Pär Johannesson, 27-Apr-99
% Toolbox: Rainflow Cycles for Switching Processes V.1.1, 27-Apr-99

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 28-Jul-2000
%   New implementation of discretization.
%   Implemented 'ddef' different methods for discretization.
% Updated by PJ 18-May-2000
%   Output 'tp' now works.
% Updated by PJ (Pär Johannesson) 25-Feb-2000
%   help
% Revised by PJ (Pär Johannesson) 12-Jan-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1999
%   Copyright (c) 1999 by Pär Johannesson, 27-Apr-99
%   Toolbox: Rainflow Cycles for Switching Processes V.1.1, 27-Apr-99

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(2,4,ni));

if ni < 3
  h=[];
end
if ni < 4
  ddef=[];
end


u = levels(param); % Discrete levels
delta = u(2)-u(1);

if isempty(h)
  h = delta;
end
if isempty(ddef)
  ddef = 1;
end

% Rainflow filter x, 
% Gives turning points with no rfc-ranges less than threshold h.

if h>0
  y = rfcfilter(x,h,0);  % Get rainflow filtered turning points
else
  y = rfcfilter(x,0,1);  % Get turning points
end

% Discretize turning points 'y'
% From TP 'y' to discrete TP 'dtp'

% Make discretization

[N,m] = size(y);
dtp = zeros(N,m);

a=param(1); b=param(2); n=param(3);
delta = (b-a)/(n-1);        % Discretization step
dtp = y;
dtp(:,m) = (y(:,m)-a)/delta + 1;

if ddef == 0
  dtp(:,m) = min(max(round(dtp(:,m)),1),n);
elseif ddef == +1
  if dtp(1,m)<dtp(2,m)  % First TP is a minimum
    dtp(1:2:N,m) = min(max(floor(dtp(1:2:N,m)),1),n-1); % Minumum
    dtp(2:2:N,m) = min(max(ceil(dtp(2:2:N,m)),2),n);    % Maximum
  else
    dtp(2:2:N,m) = min(max(floor(dtp(2:2:N,m)),1),n-1); % Minumum
    dtp(1:2:N,m) = min(max(ceil(dtp(1:2:N,m)),2),n);    % Maximum
  end
elseif ddef == -1
  if dtp(1,m)<dtp(2,m)  % First TP is a minimum
    dtp(1:2:N,m) = min(max(ceil(dtp(1:2:N,m)),2),n);    % Minumum
    dtp(2:2:N,m) = min(max(floor(dtp(2:2:N,m)),1),n-1); % Maximum
  else
    dtp(2:2:N,m) = min(max(ceil(dtp(2:2:N,m)),2),n);    % Minumum
    dtp(1:2:N,m) = min(max(floor(dtp(1:2:N,m)),1),n-1); % Maximum
  end
else
  error(['Undefined discretization definition, ddef = ' num2str(ddef)]);
end

dtp = rfcfilter(dtp,0,1);  % Get turning points

if no>=2
  uu = u';
  if m == 1
    tp = uu(dtp);
  else
    tp = [dtp(:,1) uu(dtp(:,2))];
  end
end


