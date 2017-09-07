function fout = fwaitbar(x,name,msg)
%FWAITBAR Fast display of wait bar.
%
%  CALL: H = fwaitbar(X,title,msg) 
%            fwaitbar(X,H,msg)
%
%  H     = Handle to waitbar figure.
%  X     = Fractional length of wait bar. X should be between 0 and 1.
%  title = Title string (default 'Please wait...'). 
%  msg   = Message string  (default '').
%
% FWAITBAR creates and displays a wait bar of
% fractional length X. Each subsequent call to waitbar, WAITBAR(X,H),
% extends the length of the bar to the new position X.
%
% FWAITBAR is a much speeded up version of WAITBAR (see help on WAITBAR).
%
% Example:
%  h = fwaitbar(0,[],'this may take a while');
%  for i=1:10,
%        % computation here %
%     if i==7,
%        fwaitbar(i/10,h,' soon finished')
%     else
%        fwaitbar(i/10,h)
%     end
%     pause(1)
%  end
%  close(h)
% 
% See also  waitbar

% Tested on: matlab 5.2
% History
% revised pab 03.03.2003
%  - fixed some bugs  
% revised pab 25.07.2001
% -changed help header to wafo style
% -changed old name to msg. name is now the window name. 
% - added default string to name
% - added example
% by Olof Liungman.
%   Dept. of Oceanography, Earth Sciences Centre
%   Goteborg University, Sweden
%   E-mail: olof.liungman@oce.gu.se

x = max(0,min(100*x,100)); % Make sure 0<=x<=100
if nargin<3||isempty(msg), msg = '';end
if nargin<2||isempty(name), name = 'Please wait...';end

if ischar(name)
  oldRootUnits = get(0,'Units');
  set(0,'Units','pixels');
  screenSz = get(0,'ScreenSize');
  width  = 360;
  height = 75;
  x0 = (screenSz(3)-width)/2;
  y0 = (screenSz(4)-height)/2;
  pos = [x0, y0, width, height];
  f = figure('MenuBar','none',...
     'Units','Pixels',...
     'NumberTitle','off',...
     'Pointer','watch',...
     'Color','w',...
     'Resize','on',...
     'CreateFcn','', ...
     'IntegerHandle','off',...
     'Tag','TMWWaitbar',...
     'Visible','off',...
     'Position',pos, ...
     'Name',name);
          
          %set(f,'Resize','on','Position',pos);
          %set(f,'Resize','off')
  if ispc %~strcmp(computer,'PCWIN')
    set(f,'DefaultTextFontSize',12)
    set(f,'DefaultAxesFontSize',12)
  end
  colormap([])

  ax = axes('XLim',[0 100],'YLim',[0 1],'Box','on','Position',...
       [.05 .30 .9 .30],'YTick',[],'XColor','k','YColor','k');
  
  xpatch = [0 x x 0];
  ypatch = [0 0 1 1];
  p = patch(xpatch,ypatch,'r','Edgecolor','r');%,'EraseMode','none');
  titleHandle = get(ax,'Title');
  ud = {p, titleHandle};
  set(f,'UserData',ud,...   
     'HandleVisibility','callback',...
     'Visible','on',...
     'Resize','off');
  
  set(0, 'Units', oldRootUnits);
else
   ud = get(name,'UserData');
   p           = ud{1};
   titleHandle = ud{2};

   xpatch = get(p,'xdata');
   xpatch(2:3) = x;
   set(p,'Xdata',xpatch(:))

end
if ~isempty(msg),set(titleHandle,'string',msg,'Color','k'),end

drawnow

if nargout>0,
  fout = f;
end




