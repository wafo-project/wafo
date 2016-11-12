function figpile(varargin)
%FIGPILE Pile figure windows
%
% CALL: figpile(figs)
%
%   figs = vector of integers specifying which figures to pile.
%
%   FIGPILE pile all open figure windows on top of eachother
%   with complete overlap. PILEFIGS(FIGS) can be used to specify which
%   figures that should be piled. Figures are not sorted when specified.
%
% Example:
% for ix = 1:10, figure(ix),end
% figpile             % pile all open figures
% figpile 1:3 5 7     % pile figure 1,2,3,5 and 7
%
% close all;
%  
% See also figcycle, figkeep, figmaximize, figrestore, figstack,
%          figsort, figtile

% By pab 05.03.2003
  
% based on tilefigs and stackfigs

figs = [];
if nargin>0,
  for ix=1:nargin
    if isnumeric(varargin{ix}) 
      figs = [figs,varargin{ix}];
    else
      try
        tmp = double(varargin{ix});
        if any( double('0')<= tmp & tmp<= double('9') )
          figs = [figs, eval(varargin{ix})];
        end
      end
    end
  end
end
if isempty(figs)               % If no input arguments...
   figs = findobj('Type', 'figure');    % ...find all figures.
   figs = sort(figs);                   % sort figure handles
end
if isempty(figs)
   disp('No open figures or no figures specified.');
   return
end

rootUnits = get(0, 'Units');
set(0, 'Units', 'pixels');              % Set root units.

newfigpos = get(0,'DefaultFigurePosition');

set(figs,'Position',newfigpos);

set(0, 'Units', rootUnits);  % reset root units 

for ix = figs(:).'
   figure(ix);   % raise figures
end
return


