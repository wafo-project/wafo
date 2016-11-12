function figcycle(varargin)
%FIGCYCLE Cycle through figure windows.
%
%   FIGCYCLE brings up all open figure in ascending order and pauses after
%   each figure.  When done, the figures are sorted in ascending order.
%
%   FIGCYCLE MAXIMIZE does the same thing, except that figures are maximized.
%   FIGCYCLE PAIRS2   cycle through all figures in pairs of 2.
%
% Examples:
%  for ix=1:4,figure(ix),contourf(peaks(30)),end  
%  % figcycle(1:3)          %Cycle trough figure 1 to 3
%  % figcycle 1:3           %Cycle trough figure 1 to 3 
%  % figcycle 1:3 maximize  %Cycle trough figure 1 to 3 with figs maximized   
%  % figcycle          %Cycle through all figures one at a time
%  % figtile pairs2   
%  % figcycle pairs2  % Cycle through all figures two at a time 
%
%  close all;
%
% See also figkeep, figmaximize, figrestore, figpile, figstack,
%          figsort, figtile

% Revised pab Nov 2007
%  -added option pairs
%  -added possibility to cycle backwards and escape
% Revised pab: added option 'maximize', figNumbers
%              added examples  
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-06-07 13:31:24
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://home.online.no/~pjacklam

error(nargchk(0,inf,nargin))
maximize = false;
figs = [];
nfigspercycle=1;
for ix=1:nargin
  if isnumeric(varargin{ix}) 
    figs = [figs,varargin{ix}];
  else
    try
      tmp = (varargin{ix});
      if strncmpi(tmp,'pairs',5)
        nc = numel(tmp);
        tmp = str2double(tmp(6:nc));
        if ~isempty(tmp) && tmp>0
          nfigspercycle = tmp;
        end
      elseif any( '0'<= tmp & tmp<= '9' )
        figs = [figs, eval(varargin{ix})];
      else
        maximize = true;
      end
    catch
      maximize = true;
    end
  end
end
if isempty(figs)
  % Find all figure handles, sort them and count them.
  figs = findobj('Type', 'figure');
  figs = sort(figs);
end

if isempty(figs)
   disp('No open figures or no figures specified.');
   return
end 
if maximize
  nfigspercycle = 1;
end
n    = length(figs);
%nlayers = ceil(n/nfigspercycle);

% Bring one figure up at a time.
i = 1;
while 1<=i && i<=n
  %for i = 1:nfigspercycle:n
  if maximize, %save old position and maximize window
    oldposition = get(figs(i),'position');
    figmaximize(figs(i));
  end
  for ij = 0:nfigspercycle-1
    if i+ij<=n
      figure(figs(i+ij));
    end
  end
  if i + nfigspercycle-1< n
    txt = sprintf(' %d,',figs(i + nfigspercycle:min(i+2*nfigspercycle-1,n)));
    txt(end) ='';
    fprintf('Press escape to quit, backspace to display previous figure and any other key to display figure %s\n', txt);
  end
  [xc,yc,B] = ginput(1);
  
  if maximize, %restore window position
    set(figs(i),'Position', oldposition)
  end

  if B==8, % Back space
    i = i-nfigspercycle;
  elseif B==27 % escape
    break
  else
    i = i+nfigspercycle;
  end
%    if (keydown == 0)
%      disp('Mouse button was pressed');
%    else
%      disp('Key was pressed');
%    end
end

% Sort the figures.
for i = n:-1:1
   figure(figs(i));
end

end % figcycle

