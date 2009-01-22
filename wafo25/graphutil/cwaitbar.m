function fout = cwaitbar(x,name,col)
%CWAITBAR Display compound wait bar.
%   H = CWAITBAR(X,TITLE) creates and displays wait bars of 
%   fractional lengths X and with one title string TITLE.
%   The handle to the compound waitbar figure is returned in H.
%   X values should be between 0 and 1.
%   Each subsequent call to cwaitbar, CWAITBAR([BAR X]),
%   extends the length of the bar BAR to the new position X.
%   The first bar is the topmost bar and is BAR = 1 which
%   corresponds to the outermost loop.
%   H = CWAITBAR(X,TITLE) where TITLE is a cellstring with same
%   number of titles as there are fractional lengths in X.
%   Suitable for showing the bars' corresponding loop indices.
%   H = CWAITBAR(X,TITLE,COLOR) where COLOR is the color of the
%   bars. COLOR is either a color code character (see PLOT) or
%   an RGB vector. The default color is red. COLOR can also be 
%   a cell array with same number of elements as there are bars
%   in the cwaitbar figure.
%
%   The order of the elements in vector X and cell arrays TITLE
%   and COLOR which is consistent with the bar number BAR is:
%   The first element corresponds to the first bar at the top
%   of the figure which in turn corresponds to the outermost loop.
%
%   CWAITBAR is typically used inside nested FOR loops that
%   performs lengthy computations.
%
%      Examples:
%         cwaitbar([.3 .2 .7],'Please wait...');     %single title
%
%         h = cwaitbar([0 0 0],{'i','j','k'},{[.8 .2 .8],'b','r'});
%         for i=1:5,
%            % computations %
%            for j=1:10
%               % computations %
%               for k=1:100
%                  % computations %
%                  cwaitbar([3 k/100])
%               end
%               cwaitbar([2 j/10])
%            end
%            cwaitbar([1 i/5])
%         end
%         close(h)
%
%   See also WAITBAR.

% Based on matlab's WAITBAR. See help for WAITBAR.
% Copyright (c) 2003-11-02, B. Rasmus Anthin.
% Revision 2003-11-03 - 2003-11-06.
% GPL license.

xline = [100 0 0 100 100];
yline = [0 0 1 1 0];


switch nargin
case 1   % waitbar(x)    update
   bar=x(1);
   x=max(0,min(100*x(2),100));
   f = findobj(allchild(0),'flat','Tag','CWaitbar');
   if ~isempty(f), f=f(1);end
   a=sort(get(f,'child'));                         %axes objects
   if isempty(f) || isempty(a), 
      error('Couldn''t find waitbar handles.'); 
   end
   bar=length(a)+1-bar;        %first bar is the topmost bar instead
   if length(a)<bar
      error('Bar number exceeds number of available bars.')
   end
   for i=1:length(a)
      p(i)=findobj(a(i),'type','patch');
      l(i)=findobj(a(i),'type','line');      
   end
   %rewind upper bars when they are full
%   if bar==1
%      for i=2:length(a)
%         xpatchold=get(p(i),'xdata');
%         xold=xpatchold(2);
%         if xold==100
%            set(p,'erase','normal')
%            xpatch=[0 0 0 0];
%            set(p(i),'xdata',xpatch,'erase','none')
%            set(l(i),'xdata',xline)
%         end
%      end
%   end
   
   a=a(bar);
   p=p(bar);
   l=l(bar);
   xpatchold=get(p,'xdata');
   xold=xpatchold(2);
   if xold>x                      %erase old patches (if bar is shorter than before)
      set(p,'erase','normal')
      %xold=0;
   end
   xold=0;
        %previously: (continue on old patch)
   xpatch=[xold x x xold];
   set(p,'xdata',xpatch,'erase','none')
   set(l,'xdata',xline)

case 2   % waitbar(x,name)  initialize
   x=fliplr(max(0,min(100*x,100)));

   oldRootUnits = get(0,'Units');
   set(0, 'Units', 'points');
   pos = get(0,'ScreenSize');
   pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
   
   L=length(x)*.6+.4;
   width = 360 * pointsPerPixel;
   height = 75 * pointsPerPixel * L;
   pos = [pos(3)/2-width/2 pos(4)/2-height/2 width height];

   f = figure(...
      'Units', 'points', ...
      'Position', pos, ...
      'Resize','off', ...
      'CreateFcn','', ...
      'NumberTitle','off', ...
      'IntegerHandle','off', ...
      'MenuBar', 'none', ...
      'Tag','CWaitbar');
   colormap([]);
   
   for i=1:length(x)
      h = axes('XLim',[0 100],'YLim',[0 1]);
      if ~iscell(name)
         if i==length(x), title(name);end
      else
         if length(name)~=length(x)
            error('There must be equally many titles as waitbars, or only one title.')
         end
         title(name{end+1-i})
      end
      set(h, ...
         'Box','on', ...
         'Position',[.05 .3/L*(2*i-1) .9 .2/L],...
         'XTickMode','manual',...
         'YTickMode','manual',...
         'XTick',[],...
         'YTick',[],...
         'XTickLabelMode','manual',...
         'XTickLabel',[],...
         'YTickLabelMode','manual',...
         'YTickLabel',[]);
      
      xpatch = [0 x(i) x(i) 0];
      ypatch = [0 0 1 1];
      
      patch(xpatch,ypatch,'r','edgec','r','erase','none')
      line(xline,yline,'color','k','erase','none');
      
   end
   set(f,'HandleVisibility','callback');
   set(0, 'Units', oldRootUnits);
   
case 3
   if iscell(col) && length(col)~=length(x)
      error('There must be equally many colors as waitbars, or only one color.')
   end
   f=cwaitbar(x,name);
   a=get(f,'child');
   p=findobj(a,'type','patch');
   l=findobj(a,'type','line');
   if ~iscell(col)
      set(p,'facec',col,'edgec',col)
   else
      for i=1:length(col)
         set(p(i),'facec',col{i},'edgec',col{i})
      end
   end
   set(l,'xdata',xline')
end  % case
drawnow
figure(f)

if nargout==1,
  fout = f;
end


