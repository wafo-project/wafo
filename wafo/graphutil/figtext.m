function H = figtext(x,y,T,x_unit,y_unit,x_justify,y_justify,varargin)
% FIGTEXT Places text in figure window.
%         Writes  T  at the point by  (x,y). The default coordinate
%         system used is that of the current axis. Normalized coordinates
%         are optional. Different coordinate systems for the x- and y-axis
%         are also allowed.
%
%  CALL:  H = figtext(x,y,T,x_unit,y_unit,x_justify,y_justify)
%
%        H = handle to the line of text.
%        x       = the x-coordinate of the text.
%        y       = the y-coordinate of the text.
%        T       = a string containing the text.
%        
%        x_unit,
%        y_unit  = 'data'   same coordinate system as the plotted data (default).
%                  'normalized'   use 'normalized' coordinates in [0..1].
%      x_justify = 'left'   places left adjusted text (default).
%                  'center' places centered text.
%                  'right'  places right adjusted text.
%
%      y_justify = 'top'    places top adjusted text 
%                  'middle' places centered text. (default).
%                  'bottom' places bottom adjusted text.
%
% Example: 
%   figure(1);
%   H =  figtext(0,0,'test','normalized',[],'left','top');
%
%   close all
%
%  See also  text

% Tested on: Matlab 5.3, 5.2, 5.1
% History;
% revised by pab 06.02.2000
%  added varargin
% revised by pab 11.08.99
% added y_justify

if (nargin<7)||isempty(y_justify),
  y_justify='middle';
end

if (nargin<6)||isempty(x_justify),
  x_justify='left';
end
if (nargin<4)||isempty(x_unit),
  x_unit='data';
end

if (nargin<5)||isempty(y_unit),
  y_unit=x_unit;
end
  
h=text('String',T);
set(h,'HorizontalAlignment',x_justify,'VerticalAlignment',y_justify);
unit=get(h,'Units');

set(h,'Units',x_unit)
xy=get(h,'Position');
set(h,'Position',[x xy(2)])

set(h,'Units',y_unit)
xy=get(h,'Position');
set(h,'Position',[xy(1) y])

set(h,'Units',unit);

nin=length(varargin);
if nin>0
  set(h,varargin{:});
end

if nargout>0
  H=h;
end
