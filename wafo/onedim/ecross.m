function t0 = ecross(t,f,ind,v)
%ECROSS Extracts exact level v crossings 
% 
%  CALL t0 = ecross(t,f,ind,v);
%
%  t0  = vector of  exact level v crossings.
%  t,f = vectors of arguments and functions values, respectively.
%  ind = indices to level v crossings as found by findcross.
%  v   = scalar or vector (of size(ind)) defining the level(s) to cross.
%
% ECROSS interpolates t and f linearly to find the exact level v
% crossings, i.e., the points where f(t0) = v
%
% Example
%  t = linspace(0,7*pi,250); x = sin(t);
%  ind = findcross(x,0.75);  
%  t0 = ecross(t,x,ind,0.75);
%  plot(t,x,'.',t(ind),x(ind),'r.',...
%       t, ones(size(t))*.75, t0,ones(size(t0))*0.75,'g.');
%
%  assert(ind', [10, 26, 81, 98, 152, 169, 224, 240], eps);
%  close all;
%
% See also  findcross

% Tested on: Matlab 6.0
% revised pab Feb2004  
% By pab 18.06.2001
  
error(nargchk(4,4,nargin))

t0 = t(ind)+(v-f(ind)).*(t(ind+1)-t(ind))./(f(ind+1)-f(ind));
return

