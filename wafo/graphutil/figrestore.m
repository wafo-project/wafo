function figrestore(varargin)
%FIGRESTORE  Restore figures window size and position to its default value.
%
% CALL:  figrestore(hfig)
%
%  hfig  = handle(s) of figure(s) you wish to restore (Default = 'all')
%              
%Examples:    
%  for ix = 1:5, figure(ix),end
%  figrestore('all')   %Restores all unhidden figures
%  figrestore          %same as figrestore('all')
%  figrestore(gcf)     %Restores the current figure
%  figrestore(3)       %Restores figure 3
%  figrestore([2 4])   %Restores figures 2 and 4
%  % or alternativel
%  figrestore 2 4  
%
%See also figcycle, figkeep, figmaximize, figpile, figstack,
%          figsort, figtile

%Other m-files required: none
%Subfunctions:  none
%MAT-files required:  none
%
%See also: MAXWINDOW

% revised pab 05.03.2003
% - changed default from gcf to 'all'  
% - renamed from reswindow to figrestore.
%revised pab 22.08.2002
% -updated helpheader to wafo style
%Author: Denis Gilbert, Ph.D., physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%December 2001; Last revision: 18-Feb-2002

figs = [];
for ix=1:nargin
  currArg = varargin{ix};
  if isnumeric(currArg) 
    figs = [figs,currArg];
  else
    tmp = double(currArg);
    try
      if any( double('0')<= tmp & tmp<= double('9') ) || strcmpi(currArg,'gcf')
        figs = [figs, eval(currArg)];
      end
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

rootUnits = get(0, 'Units');
set(0, 'Units', 'pixels');              % Set root units.

newfigpos = get(0,'DefaultFigurePosition');

set(figs,'units','pixels','Position', newfigpos);

set(0, 'Units', rootUnits);  % reset root units 
return
