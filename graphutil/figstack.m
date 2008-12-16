function figstack(varargin)
%FIGSTACK Stack figure windows
%
% CALL: figstack(figs)
%
%   figs = vector of integers specifying which figures to stack.
%
%   FIGSTACK stacks all open figure windows on top of eachother
%   with maximum overlap. FIGSTACK(FIGS) can be used to specify which
%   figures that should be stacked. Figures are not sorted when specified.
%
% Example:
% figstack             % stack all open figures
% figstack 1:3 5 7     % stack figure 1,2,3,5 and 7
%  
% See also: figcycle, figkeep, figmaximize, figrestore, figpile,
%          figsort, figtile

% revised pab 05.03.2003
% - adder set / reset rootUnits
% Revised pab 20.06.2001
% - updated help  header
% - added figs to input.
  

%Charles Plum                    Nichols Research Corp.
%<cplum@nichols.com>             70 Westview Street
%Tel: (781) 862-9400             Kilnbrook IV
%Fax: (781) 862-9485             Lexington, MA 02173
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


maxpos  = get (0,'screensize'); % determine terminal size in pixels

numfigs = length(figs);        % number of open figures
maxfigs = fix(maxpos(4)/20);


if (numfigs>maxfigs)            % figure limit check
        disp([' More than ' num2str(maxfigs) ' requested '])
        return
end


% tile figures by postiion 
% Location (1,1) is at bottom left corner

for iy = 1:numfigs
  figure(figs(iy))
  
  p = get(figs(iy),'Position'); % get figure position
  
  ypos = maxpos(4) - (iy-1)*20 -p(4) -70 ; % figure location (row)
  xpos = fix((iy-1)*5 + 15);     % figure location (column)
  
  set(figs(iy),'Position',[ xpos ypos p(3) p(4) ]); % move figure
end

set(0, 'Units', rootUnits);  % reset root units 
return


