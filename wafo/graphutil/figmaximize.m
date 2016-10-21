function figmaximize(varargin)
%FIGMAXIMIZE  Maximize figure(s) window size
%
% CALL:   figmaximize(hfigs, taskbar_position)
%
%  hfigs            = handle(s) of figure(s) you wish to resize (Default = 'all')
%  taskbar_position = Position of the Windows taskbar
%                     Allowed values: 'bottom',  'left'
%                     (Default = 'bottom')
%
% Note: Inputs may be given in any order.  
%  
%Examples: % Windows taskbar at bottom of screen
%  for ix = 1:5, figure(ix),end
%  figmaximize('all')   %Maximizes all unhidden figures
%  figmaximize          %same as figmaximize('all')
%  figmaximize(gcf)     %Maximizes the current figure 
%  figmaximize(3)       %Maximizes figure 3
%  figmaximize([2 4])   %Maximizes figures 2 and 4
%  figmaximize(gcf,'left')  %Windows taskbar at left of screen 
% %or alternatively
%  figmaximize 2 4  
%
% See also figcycle, figkeep, figrestore, figpile, figstack,
%          figsort, figtile


% revised pab 05.03.2003
% renamed from maxwindow to figmaximize
% changed default from gcf to 'all'  
  
%Author: Denis Gilbert, Ph.D., physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%October 1998; Last revision: 18-Feb-2002
  
taskbar_position = 'bottom'; % default

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
      else 
        switch lower(currArg(1))
          case 'b', taskbar_position = 'bottom';
          case 'l', taskbar_position = 'left';
        end
      end
    catch
    end % try-catch
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
screensize=get(0,'screensize');

switch lower(taskbar_position)
   %Other cases could easily be added for other positions of the Windows
   %taskbar or in cases where one also chooses to display the Microsoft
   %Office taskbar
   
case 'bottom'
    offset1 = 34;
    offset2 = 111;
    newPosition =[1 offset1 screensize(3) screensize(4)-offset2];
case 'left'
    offset1 = 84;
    offset2 = 78;
    newPosition = [1+offset1  1  screensize(3)-offset1  screensize(4)-offset2];
end
set(figs,'units','pixels','Position',newPosition )
set(0, 'Units', rootUnits);  % reset root units 
return