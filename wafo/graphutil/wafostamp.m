function [H,ax]=wafostamp(varargin)
%WAFOSTAMP Prints a caption "WAFO" in current figure.
%
%  CALL: [h,ax] = wafostamp(caption,flag,stamp);
%
%     h,ax  = handles to the lines of text and the axes, respectively.
%     stamp = 0, do not print WAFO stamp
%             1, print WAFO stamp         
%             2, print WAFO stamp without date (Default)
%   caption = string caption before 'made by WAFO'  (default '')
%   flag    = string following 'made by WAFO'        (default [])
%             otherwise it is customary to give one of the following flags:
%             '(ER)' - figure is Easily Reproducible ( < 10min to make)
%             '(CR)' - figure is Conditionally Reproducible ( > 10min to make)
%             '(NR)' - figure is Non Reproducible 
%
%  WAFOSTAMP creates a new invisible axes at the bottom of figure and 
%  prints the following text in it:
%  "caption"                               WAFO "date" "flag"
%  
%  NOTE: 
%        - The handles to the lines of text and the axes may also be found by 
%           h  = findobj(gcf,'tag','WAFOSTAMP','type','text');
%           ax = findobj(gcf,'tag','WAFOSTAMP','type','axes');
%  Edit wafostamp.m to change default values.
%
% Example:
%   plot(sin(0:.1:3)); wafostamp('Sinus plot','(ER)'); hold on;
%   plot(sin((0:.1:3)+pi/2)); hold off;
%
%   close all;
%
% See also  figtext

% TODO % may be further improoved by making it work like legend without the box

% Tested on matlab 5.2
% history:
% revised pab march 2007
% - removed point-and-click editing of all the text objects, because it is
%   obsolete
% revised pab 22.05.2000 minor changes
% revised pab 05.02.2000
%  -added deleteproxy
% revised pab 28.01.2000
%  - updated help header
%  - Now writes the text in a separate axes in order to make the printed
%    text independent of the scaling of the figure etc.
%  - added tag to the wafostamp axes
%  - added  point-and-click editing of all the text objects (title, xlabel, ylabel) of the current figure:
% revised pab 21.12.1999
%  - added an extra plot(sin...) in example
% revised pab 20.12.1999
%   -added varargin (caption flag), axes, date
%   -returns matlab to the same state as it was before wafostamp was called  
% adopted from watstamp by ???


[stamp,caption,flag]=stmpchk(varargin);
if ismatlab
  oldMatlabVersion = verLessThan('matlab','8.4');
 
  if ~oldMatlabVersion,
      return
   end
end
if stamp~=0
  tag = 'WAFOSTAMP';
  
  cax  = gca;
  cfig = gcf;
 
  % create new axes for caption and  delete old caption if it exists
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ax = findobj(cfig,'tag',tag);
  if ~isempty(ax), 
    delete(ax),
  end % delete old wafostamp if it exists
  ax = axes('Position',[.01 .01 .99 .05],'Visible','off','tag',tag);
 
  txtProp = {'units','normalized','tag',tag,'FontSize',10,...
    'HorizontalAlignment','left','VerticalAlignment','bottom'};
  
  
  if isempty(caption)
    h1=[];
  else
    h1 = text(0,0,caption,txtProp{:});
    %h1=figtext(0,0,caption,'normalized','normalized','left','bottom');
    %set(h1,'FontSize',10,'Tag',tag)
  end
  if 1,
    
    bannerFile = fullfile(waforoot,'data','wafologoNewWithoutBorder.png');
    %setbanner(handles,bannerFile)
    [banner,MAP,ALPHA] = imread( bannerFile); % Read the image file banner.jpg
    %info           = imfinfo(bannerFile); % Determine the size of the image file
  
    ax(2) = axes('Position',[.70 .01 .2 .05],'Visible','off','tag',tag);
    axis equal
 
    axes(ax(2));
  
    image(banner)
    set(ax(2), 'Visible', 'off'); %, 'Position', [50 50 info.Width info.Height]);
    if stamp>1
      stamptxt =[flag,' '];
    else
      stamptxt =[date,' ',flag,' '];
    end
    
  else
    % Old call
    if stamp>1
      stamptxt =['WAFO ',flag,' '];
    else
      stamptxt =['WAFO ' date,' ',flag,' '];
    end
  end
   txtProp1 = {'units','normalized','tag',tag,'FontSize',10,...
     'HorizontalAlignment','right','VerticalAlignment','bottom','FontAngle','Italic'};
  axes(ax(1));
   h2 = text(1,0,stamptxt,txtProp1{:});
%    h2=figtext(1,0,stamptxt,'normalized','normalized','right','bottom');
%   % old call
%   %h=figtext(0,-0.1,'made by WAFO','normalized','normalized','left');
%   set(h2,'FontSize',10,'FontAngle','Italic','Tag',tag)
  
  
  
  % create DeleteProxy object (an invisible text object in
  % the first axes) so that the other axes will be deleted
  % properly.
  %'tag','wafostmptxt',...
  mkdeleteproxy(cax,ax,'WAFOSTAMPDeleteProxy');


  if nargout>0
    H=[h1 h2];
  end

  % reset to the old state
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %set(cfig,'currentaxes',cax)  % reset to old axes
  axes(cax);
  
  
end


function [stamp,caption,flag]=stmpchk(P)
% gives the correct input values  
N = length(P);

stamp   = 2; %EDIT HERE to change default value
caption = '';
flag    = '';
if N>0
  isNum = cellfun(@isnumeric,P);

  k = find(isNum);
  if any(k) && ~isempty(P{k})
    stamp = P{k};
  end
  k1 = find(~isNum);
  if any(k1)
    caption = P{k1(1)};
    if length(k1)>1
      flag = P{k1(2)};
    end
  end
end

  
