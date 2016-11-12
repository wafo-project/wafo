function [H] = cltext(clevels,PL,N)
%CLTEXT Places contour level text in the current window
%
% CALL: h = cltext(levels,pl,N)
%            
%      h       = handles to the text objects.
%      levels  = vector of contour levels or the corresponding percent which the
%                contour line encloses
%      pl      = 0 if levels are the actual contour levels (default)
%                1 if levels are the corresponding percent which the
%                  contour line encloses
%      N       = maximum N digits of precision (default 4)
%
% CLTEXT creates text objects in the current figure and prints 
%        "Level curves at:"        if pl==0 and
%        "Level curves enclosing:" otherwise
%  and the contour levels or percent.
%
%  NOTE: 
% -The handles to the lines of text may also be found by 
%        h  = findobj(gcf,'tag','CLTEXT','type','text');
%        h  = findobj(gca,'tag','CLTEXT','type','text');
% -To make the text objects follow the data in the axes set the units 
%  for the text objects 'data' by    
%        set(h,'unit','data')
%
% Examples:
%  z  = peaks;
%  cl = max(z(:))-range(z(:))*(.1:.1:.9);
%  contour(z,cl); cltext(cl);
%
%  data = rndray(1,2000,2); 
%  f = kdebin(data,{'kernel','epan','L2',.5,'inc',128});
%  contour(f.x{:},f.f,f.cl); cltext(f.pl,1);
%
%  close all;
%
% See also  pdfplot



% History
% Revised pab March 2007
% - The printed text is now truly independent of the scaling to the rest of
%   the figure.
% revised pab Feb2004
%  made it work with zoom to some extent  
% revised pab Jan2004
% changed todo comments: (+) -> TODO  
% revised pab 22.05.2000 minor changes
% revised pab 28.01.2000
% - prints the text in a separate axes in order to make the printed
%    text independent of the scaling of the rest of the figure etc.
% adopted from pdfplot

% Copyright (C) 2000  Per A. Brodtkorb
% 
%  This file, CLTEXT.M, is part of WAFO.
% 
%     CLTEXT is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     CLTEXT is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    


% TODO % Make it work like legend does (but without the box): include position options etc...


if nargin<1||isempty(clevels)
  return % nothing to do
end
if nargin<2||isempty(PL)
  PL=0;
elseif numel(PL)>1
  error('pl must bea scalar: 0 or 1')
end
if nargin<3 || isempty(N)
  N = 4;
end

tag = 'CLTEXT';
cax    = gca; % get current axes


% delete cltext object if it exists 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_CLTXT = findobj(cax,'tag',tag); 
if ~isempty(H_CLTXT)
  try
    delete(H_CLTXT)
  catch
    warning('WAFO:CLTEXT','Tried to delete a non-existing CL-text')
  end
end


% New call that does not work properly
%oldUnit = get(cax,'FontUnits');
%set(cax,'FontUnits','normalized');
%charHeight = get(cax,'FontSize');
%set(cax,'FontUnits',oldUnit);

charHeight = 1/33;
delta_y     = charHeight;

  % [textstart_x,textstart_y,1,1]
if (PL==1),
  titletxt = 'Level curves enclosing:';
else
  titletxt = 'Level curves at:';
end




format = ['%0.' int2str(N),'g\n'];


if 1
  cltxt = sprintf(format,clevels);
else
  cltxt = num2str(clevels(:),N);
  
  %removing spaces in front of each line
  indx = find(isspace(cltxt(:,1)));
  space = ' ';
  for ix=indx(:).'
    ik = find(~isspace(cltxt(ix,:)),1);
    cltxt(ix,:) = [cltxt(ix,ik:end) space(:,ones(1,ik-1))];
  end
  %cltxt = cellstr(cltxt);
end


textstart_x = 0.036;
textstart_y = 0.94;
textstart_z = 0;

position = {textstart_x,textstart_y,textstart_z};

if isoctave,
   warning ('off', 'Octave:abbreviated-property-match')
end

titleProp = {'unit','normalized','tag',tag,...
  'HorizontalAlignment','left',...
  'VerticalAlignment','middle',...
  'FontWeight','bold'};
ha(1) = text(position{:},titletxt,titleProp{:});

%textstart_y = textstart_y - delta_y;
position{2} = position{2} - delta_y;

txtProp = {'unit','normalized','tag',tag,...
  'HorizontalAlignment','left','VerticalAlignment','top'};
ha(2) = text(position{:},cltxt,txtProp{:});

% create DeleteProxy object, an invisible text object in
% the current axes, so that the other text objects will be deleted
% properly.
mkdeleteproxy(cax,ha,'CLTEXTDeleteProxy')

%set(ha,'unit','data')

if nargout>0
  H = ha;
end


return %% main



