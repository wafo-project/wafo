function [y] = rfcfilter(x,h,def)  
% RFCFILTER Rainflow filter a signal.
%
% CALL:  [y] = rfcfilter(x,h,def);
% 
% Input:
%   x   = Signal.   [nx1] OR [nx2]
%   h   = Threshold for rainflow filter.
%   def = 0: removes cycles with range < h. (default)
%         1: removes cycles with range <= h.
% Output:
%   y   = Rainflow filtered signal.
%
% Examples:  % 1. Filtered signal y is the turning points of x.
%      x = load('sea.dat');
%      y = rfcfilter(x,0,1);
%            % 2. This removes all rainflow cycles with range less than 0.5.
%      y = rfcfilter(x,0.5);

% Tested on Matlab 6.0
%
% History: 
% Corrected final check, added lines 161-165, by GL, 25-Jan-2018
% Updated by PJ 18-May-2000
%   Help text.
% Revised by PJ 12-Jan-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1999

  
% Check input and output

ni = nargin;
no = nargout;
%error(nargchk(2,3,ni));
narginchk(2,3)
if ni < 3
  def=[];
end

% Set default values

if isempty(def)
  def = 0;
end

% Initiation

[n,m] = size(x);
y = zeros(n,m);

y(1,:) = x(1,:);
j = 1;
if m == 1
  t0 = 0;
  y0 = y(1);
else
  t0 = y(1,1);
  y0 = y(1,2);
end

z0 = 0;

% The rainflow filter

for i = 2:n
  fpi = y0+h; %fp(y0,h);
  fmi = y0-h; %fm(y0,h);
  if m == 1
    ti = 0;
    xi = x(i);
  else
    ti = x(i,1);
    xi = x(i,2);
  end
  if z0 == 0
    if def == 0
      test1 = (xi <= fmi);
      test2 = (xi >= fpi);
    else % def == 1
      test1 = (xi < fmi);
      test2 = (xi > fpi);
    end
    if test1
      z1 = -1;
    elseif test2
      z1 = +1;
    else
      z1 = 0;
    end
    if z1 == 0
      t1 = t0;
      y1 = y0;
    else
      t1 = ti;
      y1 = xi;
    end
  else
    
    % Which definition?
    if def == 0
      test1 = (((z0==+1) & (xi<=fmi))  | ((z0==-1) & (xi<fpi)));
      test2 = (((z0==+1) & (xi>fmi)) | ((z0==-1) & (xi>=fpi)));
    else % def == 1
      test1 = (((z0==+1) & (xi<fmi))  | ((z0==-1) & (xi<=fpi)));
      test2 = (((z0==+1) & (xi>=fmi)) | ((z0==-1) & (xi>fpi)));
    end
    
    % Update z1
    if     test1
      z1 = -1;
    elseif test2
      z1 = +1;
    else
      warning(['Something wrong, i=' num2str(i)]);
    end
    
    % Update y1
    if     z1 ~= z0
      t1 = ti;
      y1 = xi;
    elseif z1 == -1
      % y1 = min([y0 xi]);
      if y0<xi
	t1 = t0;
	y1 = y0;
      else
	t1 = ti;
	y1 = xi;
      end
    elseif z1 == +1
      % y1 = max([y0 xi]);
      if y0 > xi
	t1 = t0;
	y1 = y0;
      else
	t1 = ti;
	y1 = xi;
      end
    end
    
  end
  
  % Update y if y0 is a turning point
  if abs(z0-z1) == 2
    j=j+1;
    if m ==1
      y(j) = y0;
    else
      y(j,:) = [t0 y0];
    end
  end
  
  % Update t0, y0, z0
  t0=t1;
  y0=y1;
  z0=z1;
end

% Update y if last y0 is greater than (or equal) threshold
if m == 1
    yj = y(j);
else
    yj = y(j,2);
end
if def == 0
  test = (abs(y0-yj) > h)
else % def == 1
  test = (abs(y0-yj) >= h);
end
if test
  j=j+1;
  if m ==1
    y(j) = y0;
  else
    y(j,:) = [t0 y0];
  end
end

% Truncate y
y=y(1:j,:);

function uu = fp(u,h)
  uu = u+h;
  
function uu = fm(u,h)
  uu = u-h;
