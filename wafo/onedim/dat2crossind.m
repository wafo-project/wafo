function [ind , Nc]= dat2crossind(x,v,wdef,nowarning)
%DAT2CROSSIND Finds indices to level v down and/or upcrossings from data
%
% CALL:  [ind, Nc]= dat2crossind(x,v,wdef/cdef,nowarning);
%
%  ind  = indices to the level v crossings of the original sequence x
%  Nc   = number of crossings (i.e.length of ind) 
%  x    = the surface elevation data
%  v    = the reference level (default  v = mean of  x)
%  wdef = defines the type of wave. Possible options are
%        'dw', 'uw', 'cw', 'tw' or 'none'. (Default 'none').
%        If wdef='none' all crossings will be returned,
%        otherwise only the crossings which defines a 
%        wave according to the wave definition will be returned.
%  cdef = defines the type crossings returned. Possible options are
%        'd' 'u' or 'all'. (Default 'all').
%        If def='d' all down-crossings will be returned.
%        Similarly if def='u' only the up-crossings will be returned
%        otherwise 'all' the crossings will be returned.
%  nowarning = true suppresses warning for no crossings (default = false)
%
% Example: 
%   t = linspace(0,7*pi,250); 
%   x = sin(t);
%   [ind, Nc] = dat2crossind(x,0.75,'u');
%   plot(t,x,'.',t(ind),x(ind),'o');
%
%   assert(ind', [ 10,81,152,224], eps);
%   close all;
% See also  findcross, wavedef, crossdef

%Tested on: Matlab 8.6, 8.1, 6.0, 5.3, 5.2, 5.1
% History:
% added option to suppress warning fo "no crossings found" GL 03-04-2015
% updated to Matlab 2012+, GL 31.10.2014
% suppressed empty warning GL 2006-12-28
% revised pab Feb2004
% revised by pab 12.06.2001
%  -added check on ind returned from findcross
% Revised by jr 02.04.2001
% - Added example, updated help 
% By Per A. Brodtkorb 07.07.1998,  27.07.1998.  

%error(nargchk(1,4,nargin))
narginchk(1,4)
xn=x;
if nargin<4,
    nowarning=false;
end
[n, m]= size(xn);
if n<m
 b=m;m=n;n=b; 
 xn=xn';
end

if n<2, 
  error('The vector must have more than 2 elements!');
end

%istime=1;

switch m
 case 1, %istime=0;
 case 2, xn= xn(:,2);% dimension OK!
 otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ');
end

if ((nargin<3) || isempty(wdef)),
  wdef='none';
end

if ((nargin<2) || isempty(v)),
  v = mean(xn);
  disp(['   The level v is set to: ', num2str(v)]);
end


% find level v down-crossings and/or up-crossings
% according to wdef or cdef

ind = findcross(xn,v); % faster than find

if isempty(ind), %added pab 12.06.2001
   Nc = 0; 
   if nowarning == false
       txt = sprintf('No level v = %0.5g crossings found in x',v);
       warning(txt);
   end
   return,
end


switch wdef   % switch wdef/cdef    
  case 'd', %downcrossings only
    if xn(ind(1)+1)>v,
      ind =ind(2:2:end);
    else
      ind =ind(1:2:end);
    end
    
 case 'u',%upcrossings  only
   if xn(ind(1)+1)<v,
      ind =ind(2:2:end);
    else
      ind =ind(1:2:end);
    end
    
  case {'dw','uw'},
    % make sure that the first is a level v down-crossing if wdef == 'dw'
    % or make sure that the first is a level v up-crossing if wdef == 'uw'

    if xor(((xn(ind(1))>xn(ind(1)+1))),strcmp(wdef,'dw')),
      ind(1)=[];
    end
    Nc=length(ind); % number of level v crossings
    % make sure the number of troughs and crests are according to the
    % wavedef, i.e., make sure length(ind) is odd
    if ~(mod(Nc,2)), % if Nc is even do
      ind(end)=[];
    end
  case {'tw','cw'},
    % make sure that the first is a level v down-crossing if wdef == 'tw'
    % or make sure that the first is a level v up-crossing if wdef == 'cw'

    if xor(((xn(ind(1))>xn(ind(1)+1))),strcmp(wdef,'tw')),
      ind(1)=[];
    end
    Nc=length(ind); % number of level v crossings
    % make sure the number of troughs and crests are according to the
    % wavedef, i.e., make sure length(ind) is even
    if (mod(Nc,2)), % if Nc is odd do
      ind(end)=[];
    end
  case {'du','all','none'},
    % do nothing
 otherwise,  error('Unknown wave/crossing definition!');
end
if nargout>1,
  Nc=length(ind); % number of level v crossings
end

return



