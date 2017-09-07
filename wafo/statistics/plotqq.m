function varargout=plotqq(x,y,ps,method)
%PLOTQQ   Plot empirical quantile of X vs empirical quantile of Y
%
% CALL:  h = plotqq(x,y,ps)
%
%  h   = handle to the plotted figure
%  x,y = data 
%  ps  = plot symbol (default '.')
%
% If two distributions are the same (or possibly linearly 
% transformed) the points should form an approximately straight 
% line.
%
% Example:
%   R1=rndgumb(1,0,1,100);
%   R2=rndgumb(2,2,1,100);
%   plotqq(R1,R2)
%
% See also  percentile

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Adapted to  cssmooth  by GL Feb 2011
% tested on: matlab 5.3, 7, 8, 9
% History:
% revised pab March 2007
% - added LS-line to qqplot using smooth
% - handle to plotted figure is only returned if nargout>0
% revised pab 24.10.2000
% added nargchk
% updated header to comform to wafo style

%error(nargchk(2,4,nargin))
narginchk(2,4)

if nargin<3 || isempty(ps), 
  ps = '.'; 
end
ps2 = 'r-.';
if nargin<4 || isempty(method), 
  method = 1; 
end
x = sort(x);
y = sort(y);
nx = length(x);
ny = length(y);
n = min(nx,ny);
ne = max(floor(n/10),0);
ix = ne+1:n-ne;


if nx < ny
   yi = percentile(y, ((1:n)-0.5)/n, method);
   h =plot(x, yi , ps,x,lsline(x(ix),yi(ix),x),ps2);
   
elseif ny < nx
   xi = percentile(x, ((1:n)-0.5)/n, method);
   h = plot(xi, y, ps,x,lsline(xi(ix),y(ix),xi),ps2);
else
   h=plot(x,y,ps,x,lsline(x(ix),y(ix),x),ps2);
end
xlabel('Data X')
ylabel('Data Y')
title('Quantile Plot')
wafostamp;

if nargout>0
 varargout{1:nargout} = h;
end

function  yi = lsline(x,y,xi)
% LSLINE Least square line

[x,ix] = unique(x(:)); % Remove ties.

% ...and small stepsizes

dx = diff(x);
dx = [dx(1);dx];
dxsmall = percentile(dx,0.25);
iy = find(dx>=dxsmall);
ix = ix(iy);
yi = cssmooth(x(iy),y(ix),0,xi);
