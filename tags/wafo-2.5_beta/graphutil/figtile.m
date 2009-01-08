function figtile(varargin)
%FIGTILE  Tile figure windows.
%
% CALL:  figtile(figs)
%
%   figs = vector of integers specifying which figures to tile. 
%
%   FIGTILE places all open figure windows around on the screen with no
%   overlap. FIGTILE(FIGS) can be used to specify which figures that
%   should be tiled. Figures are not sorted when specified.
%
% Example: 
% for ix=1:10,figure(ix);end
% figtile             % tile all open figures
% figtile 1:3 5 7     % tile figure 1,2,3,5 and 7
% figtile pairs2 1:10 % tile figure 1 to 10 two at a time
% figtile pairs3 1:10 % tile figure 1 to 10 three at a time     
%  
% See also figcycle, figkeep, figmaximize, figrestore, figpile, figstack, figsort


%
% revised pab Nov 2007
%  -added option pairs
% revised pab 05.03.2003
% - added reset rootUnits  
%  
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-06-07 13:28:56
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://home.online.no/~pjacklam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the handles to the figures to process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figs = [];
nfigspertile = [];
for ix=1:nargin
  currArg = varargin{ix};
  if isnumeric(currArg) 
    figs = [figs,currArg];
  else
    tmp = double(currArg);
    try
      if strncmpi(currArg,'pairs',5)
        tmp = str2double(currArg(6:end));
        if ~isempty(tmp) && tmp>0
          nfigspertile = tmp;
        end
      elseif any( double('0')<= tmp & tmp<= double('9') ) || strcmpi(currArg,'gcf')
        figs = [figs, eval(currArg)];        
      end
    catch
      % 
    end
  end
end
if isempty(figs)
  % Find all figure handles, sort them and count them.
  %figs = get(0,'children');
  figs = findobj('Type', 'figure');
  figs = sort(figs);
end

if isempty(figs)
   disp('No open figures or no figures specified.');
   return
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set miscellaneous parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfigs = length(figs);           % Number of figures.
if isempty(nfigspertile)
  nfigspertile = nfigs;
end
nlayers = ceil(nfigs/nfigspertile);

nh = ceil(sqrt(nfigspertile));         % Number of figures horisontally.
nv = ceil(nfigspertile/nh);            % Number of figures vertically.


nh = max(nh, 2);
nv = max(nv, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the screen size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rootUnits = get(0, 'Units');

set(0, 'Units', 'pixels');              % Set root units.
scrdim = get(0, 'ScreenSize');          % Get screen size.
scrwid = scrdim(3);                     % Screen width.
scrhgt = scrdim(4);                     % Screen height.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The elements in the vector specifying the position.
% 1 - Window left position
% 2 - Window bottom position
% 3 - Window horizontal size
% 4 - Window vertical size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hspc   = 50;            % Horisontal space.
topspc = 75;            % Space above top figure.
medspc = 75;            % Space between figures.
botspc = 40;            % Space below bottom figure.

figwid = ( scrwid - (nh+1)*hspc )/nh;
fighgt = ( scrhgt - (topspc+botspc) - (nv-1)*medspc )/nv;

figwid = round(figwid);
fighgt = round(fighgt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put the figures where they belong.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 0;
for ix = 1:nlayers
  for row = 1:nv
    for col = 1:nh
      %idx = (row-1)*nh + col; 
      if  (row-1)*nh + col <= nfigspertile      % idx <= nfigs
        idx = idx+1;
        if idx<=nfigs
        figlft = col*hspc + (col-1)*figwid;
        figbot = scrhgt - topspc - row*fighgt - (row-1)*medspc;
        figpos = [ figlft figbot figwid fighgt ];    % Figure position.
        fighnd = figs(idx);                          % Figure handle.
        units = get(fighnd, 'Units');                % Get units.
        set(fighnd, 'Units', 'pixels');              % Set new units.
        set(fighnd, 'Position', figpos);             % Set position.
        set(fighnd, 'Units', units);                 % reset units
        figure(fighnd);                              % Raise figure.
        end
      end
    end
  end
end
set(0, 'Units', rootUnits); 


