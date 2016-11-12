function [h00,h1]=fcolorbar(D,L,label,W,hand) 
%FCOLORBAR  Display colorbar with discrete color axis for filled contour plot
%
% CALL:  [h0,h1] = fcolorbar(C,L,label,W,hand);
%        [h0,h1] = fcolorbar(V,L,label,W,hand);
%
% V     = The contour-level-specification vector used when making the 
%         CONTOURF-plot (the most robust method), 
% C     = The contour-matrix as described in CONTOURC
% L     = String giving location of colorbar:
%         't'op
%         'b'ottom
%         'r'ight   (default)
%         'l'eft
%         'o'utside (when rescaling of graph is undesireable)
%         'a'lone   (in a subplot when it refers to several
%          graphs in a subplot-figure)
% label = string for the color-axis-label               (default= empty)
% W     = Relative width (ratio of graph width/height)  (default= 1/20)
% hand  = handles to axes which the colorbar refers to. Use this to
%	  ensure the right colorspan on all graphs and their colorbar.
%
% [h0,h1] are handles to the graph and to the colorbar axes, respectively.
%
% When contours (V) are given, the colorlimit ('Clim') property of
% both the plot and the colorbar is locked to the lower and upper
% contour level.  This should solve the problem of how to keep the same
% color-to-data relationship from plot to plot. 
%
% As a general rule, place the colorbar on a side of the graph that has
% no ticklabels or title (right is this function's default).
%
% NOTE : CONTOURF scales CAXIS from the lowest input contourlevel to the
%         highest given contourlevel or max(max(data)) whichever is
%         greatest!!! :-(
%         FCOLORBAR locks the coloraxis of both plot and colorbar between
%         the first and last contourlevel, to ensure consistency. When
%         several axes 'share' same colorbar, the 'clim' property of all
%         axes should be set equal to the colorbar, by input of their
%         handles in 'hand'! 
%
% Examples: % Difference from the ordinary colorbar-function:
%
% [x,y,z]=peaks; v=-10:2:8; 
% figure(1);  [cs,h]=contourf(x,y,z,v); clabel(cs,h); colorbar
% figure(2);  [cs,h]=contourf(x,y,z,v); clabel(cs,h); fcolorbar(v);
%           % And not using contourspecification:
% figure(3);  [cs,h]=contourf(x,y,z); clabel(cs,h); fcolorbar(cs);
%
%           % Not equidistant contours.
% figure(4);  v=[-8 -4 -2 -1 0 1 2 4 8]; 
% [cs,h]=contourf(x,y,z,v); clabel(cs,h); fcolorbar(v);
%
% close all;
%
% See also  contourf, clevels


% ### Updates: ###
% revised pab 31.12.2002
% -possible to give a single contour-level specification
% -added deleteProxy object
% -renamed from ecolorbar to fcolorbar
% -moved some code into separate functions
% revised pab 20.08.2001
% - streamlined some code.
% revised pab 03.07.2001
% - changed help header to WAFO style.
% 00.10.17: Added options for "outside"- and "alone"-positioning 
% 99.11.19: I think I found the solution for the color-span problem.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modified from my_colorbar by Tore Furevik, Geophysical Institute,
% University of Bergen. 
% ECOLORBAR by Jan Even Nilsen:
%Time-stamp:<Last updated on 01/01/30 at 12:04:59 by even@gfi.uib.no>
%File:<d:/home/matlab/ecolorbar.m>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(1,5,nargin));

% Catch fcolorbar('delete') special case -- must be called by the deleteFcn.
if nargin==1 && strcmp(D,'delete'),
  ax = gcbo;
  if strcmp(get(ax,'tag'),'ColorbarDeleteProxy')
     cbo = ax;
     ax = get(cbo,'userdata');
     if ishandle(ax)
        ud = get(ax,'userdata');
        
        % Do a sanity check before deleting colorbar
        if isfield(ud,'DeleteProxy') && isequal(ud.DeleteProxy,cbo) 
           try
              delete(ax)
           end
           if   ishandle(ud.DeleteProxy)
              try
                 delete(ud.DeleteProxy)
              end
           end
           
           if isfield(ud,'PlotHandle') && ishandle(ud.PlotHandle) && ...
                 isfield(ud,'originalPosition') && ~isempty(ud.originalPosition)
              units = get(ud.PlotHandle,'units');
              set(ud.PlotHandle,'units','normalized');
              set(ud.PlotHandle,'position',ud.originalPosition);
              set(ud.PlotHandle,'units',units);
              
              updateLegend(ud.PlotHandle); 
          
           end
        end
     end
  elseif isempty(ax),
        %fcolorbar delete call
     ax = findobj(gcf,'type','axes','tag','Colorbar');
     ud = get(ax,'userdata');
     if isfield(ud,'PlotHandle') && ishandle(ud.PlotHandle) && ...
           isfield(ud,'originalPosition') && ~isempty(ud.originalPosition)
        units = get(ud.PlotHandle,'units');
        set(ud.PlotHandle,'units','normalized');
        set(ud.PlotHandle,'position',ud.originalPosition);
        set(ud.PlotHandle,'units',units);
        
        updateLegend(ud.PlotHandle); 
     end
     if isfield(ud,'DeleteProxy') && ishandle(ud.DeleteProxy)
        try
           delete(ud.DeleteProxy)
        end
     end
    
  end
  return
end

if nargin < 5 || isempty(hand),  hand = [];end
if nargin < 4 || isempty(W),     W    = 1/20;end
if nargin < 3 || isempty(label) || ~ischar(label),  label='';end
if nargin < 2 || isempty(L),     L    = 'r';end  

outside=0; alone=0;
L = lower(L);
switch L(1)
  case 'o', outside = 1; L(1) = [];
  case 'a', alone   = 1; L(1) = [];
end
if isempty(L), L = 'r';end


% What D is input...
Dsiz = size(D);
if prod(Dsiz)==max(Dsiz) && all(Dsiz>0) % D is a vector and therefore a contourspecification
  V = unique(D);
elseif Dsiz(1)==2,
  V = clevels(D);     % contour matrix given
else
  error('First input must be a vector of contour level specifications or a contour-matrix!');
end             

%moved tick and limits here to streamline the code. pab 20.08.2001
%tick    = maketick(V);
tick    = V;
if (length(V)>1)
   mv      = mean(diff(V(:)))/2; % average stepsize divided by 2
   limits = [min(V)-mv max(V)+mv];
else
   mv = max(abs(V),1);
   limits = V + [-mv  mv];
end
logScale = 0;
if logScale,
   if min(V)<=eps
      error('log scale impossible')
   end
   if limits(1)<0
      limits(1) = eps;
   end
end
V       = [limits(1) ;V(:); limits(2)];

[x,VC]  = meshgrid(1:2,V);  % using the given contour-range 

forceEqualSpacing = 1;
if forceEqualSpacing,
   [x,y] = meshgrid(1:2, 1:length(V));
else
   y = VC;
end

                         
% h0 = handle to the graph's axes ; 
% h1 = handle to the colorbar's axes
h0 = gca; 


if alone, set(h0,'visible','off'); end

h1 = findColorbarAxes(h0);

if isempty(h1), % no colorbar present
   
  %legend('RestoreSize',h0); %restore axes to pre-legend size
  
  % Make colorbarAxes and update graph axes position
  units = get(h0,'units'); 
  set(h0,'units','normalized')
  originalPosition  = get(h0,'position');
  [cAxesPosition,gAxesPosition] = getNewAxesPositions(h0,originalPosition,L,outside,W);   
  set(h0,'Position',gAxesPosition);
  set(h0,'units',units)
  h1 = axes('Position', cAxesPosition);
  
  %Create a userData object saving important info about the colorbar, graph axes and deleteProxy 
  ud.originalPosition = originalPosition;
  ud.PlotHandle       = h0;
  ud.Location         = L(1);
  ud.outside          = outside;
  
  % create DeleteProxy object (an invisible text object in
  % the colorbar axes) so that the colorbar will be deleted
  % properly.
  deleteProxyProperties = {'parent',h0,'visible','off','tag','ColorbarDeleteProxy',...
        'handlevisibility','off','deletefcn','fcolorbar(''delete'')'};%'eval(''delete(get(gcbo,''''userdata''''))'','''')');%
  ud.DeleteProxy = text(deleteProxyProperties{:});
 
  set(ud.DeleteProxy,'userdata',h1)
  
else  % colorbar axes already exist
   ud = get(h1,'UserData');
   % Make sure contourf deletefcn doesn't trigger a fcolorbar('delete')
   % for colorbar update
   %set(get(h1,'children'),'deletefcn','')
    
   sameColorbarLocation = strcmpi(ud.Location,L(1));
   if ~sameColorbarLocation || ud.outside ~=outside,
      %update userData object
      ud.Location = L(1);
      ud.outside  = outside;
      
      % update colorbar Axes and graph axes positions
      
      %set(h0,'Position',ud.originalPosition);
      [cAxesPosition,gAxesPosition] = getNewAxesPositions(h0,ud.originalPosition,L,outside,W);
      units = get(h0,'units'); 
      set(h0,'units','normalized')
      set(h0,'Position',gAxesPosition);
      set(h0,'units',units)
      set(h1,'Position',cAxesPosition);
      
   end
   axes(h1); 
end


 % Plot the actual colorbar
if any(L(1)=='tb') 
   %[x,y] = deal(y,x); % swap x and y
   contourf(y,x,VC,V);
   if (L(1)=='t'), set(h1,'XAxisLocation','top');  end
   set(h1,'XTick',tick,'Xlim',limits,'YTick',[]); 
   set(h1,'XScale','log')
    if forceEqualSpacing
      set(h1,'YTick',[],'XTick',2:length(tick)+1,'Xlim',[1 length(V)]);
      set(h1,'Xticklabel',num2str(tick(:)));
   end
   xlabel(label)
else
   contourf(x,y,VC,V); 
   if (L(1)=='r'), set(h1,'YAxisLocation','right'); end

   set(h1,'XTick',[],'YTick',tick,'Ylim',limits);
   if logScale,
      set(h1,'YScale','log')
   end
   if forceEqualSpacing
      set(h1,'XTick',[],'YTick',2:length(tick)+1,'Ylim',[1 length(V)]);
      set(h1,'Yticklabel',num2str(tick(:)));
   end
   ylabel(label)
end


ud.PlotHandle = h0;
set(h1,'userdata',ud,'Tag','Colorbar');
set(gcf,'CurrentAxes',h0); % set current axis pointer back on the graph

%set(ud.DeleteProxy,'deletefcn','fcolorbar(''delete'')');

updateLegend(h0);

% lock the coloraxis onto the colorbar in both plot and bar:
set([h0 h1 hand(:)'],'clim',limits);

if nargout>0,  h00 = h0; end % pab 20.08.2001

return


%-------------------------------------------------------------
function tick=maketick(lev)
%MAKETICK Make the number of ticklabels be less than 10
%
% CALL: tick = maketick(levels)
tick=lev;
while length(tick)>10
  tick=tick(1:2:length(tick));
end

% no tick on the edges of colorbar
%if max(tick)==max(lev)
%  tick=tick(1:length(tick)-1);
%end
return




function cax = findColorbarAxes(dataAxes)
%FINDCOLORBARAXES Return the handle to the colorbar of dataAxes if it exists.
%

cfig = get(dataAxes,'parent'); % get figure handle

% Search for existing colorbar if any
ch = findobj(cfig,'type','axes','tag','Colorbar');
cax = [];
for i=1:length(ch),
   ud = get(ch(i),'UserData'); 
   d = ud.PlotHandle;
   if numel(d)==1 && isequal(d,dataAxes), 
      cax = ch(i);
      break; 
   end
end
return


function [cAxesPosition,gAxesPosition] = getNewAxesPositions(h0,a,L,outside,W)
%GETNEWAXESPOSITIONS Return the new normalized axes positions for the colorbar and the graph, respectively.
%
% CALL: [cAxesPosition,gAxesPosition] = getNewAxesPositions(gax,L,outside,W); 
%
% cAxesPosition
% ,gAxesPosition     = new normalized axes positions for the colorbar and the graph, respectively.
%  gax               = handle to the graph axes.
%  gAxesOrigPosition = original axes position for the graph (normalized)
%  L                 = String giving location of colorbar:
%                      't'op
%                      'b'ottom
%                      'r'ight
%                      'l'eft
% outside            = 0 when rescaling of graph is wanted
%                      1 when rescaling of graph is undesireable
% W                  = Relative width (ratio of graph width/height)
%


%a = get(h0,'Position');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch L(1)
case 't'  
   W=W*a(4);
   % placement of the colorbar:
   if findstr(get(h0,'XAxisLocation'),'top'),os=2*W; else os=W; end 
   if outside
      cAxesPosition = [a(1) a(2)+a(4)+os  a(3) W];
      gAxesPosition = a;
   else
      cAxesPosition = [a(1) a(2)+a(4)-W  a(3) W];
      gAxesPosition = [a(1) a(2) a(3) a(4)-W-os];
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'b'     
   W=W*a(4);
   % placement of the colorbar:
   if findstr(get(h0,'XAxisLocation'),'bottom'),os=W*2; else os=W; end 
   if outside
      cAxesPosition = [a(1) a(2)-os-W a(3) W];
      gAxesPosition = a;
   else
      cAxesPosition = [a(1) a(2) a(3) W];
      gAxesPosition = [a(1) a(2)+os+W a(3) a(4)-os-W];
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'r'  %default
 W=W*a(3);
 % placement of the colorbar:
 if findstr(get(h0,'YAxisLocation'),'right'),os=W*2; else os=W; end 
 if outside
    cAxesPosition = [a(1)+a(3)+os a(2) W a(4)];
    gAxesPosition = a;
 else
    cAxesPosition = [a(1)+a(3)-W a(2) W a(4)];
    gAxesPosition = [a(1) a(2) a(3)-os-W a(4)];
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'l'   
   W=W*a(3);
   % placement of the colorbar:
   if findstr(get(h0,'YAxisLocation'),'left'),os=W*2; else os=W; end 
   if outside
      cAxesPosition = [a(1) a(2)-os-W W a(4)];
      gAxesPosition = a;
   else
      cAxesPosition = [a(1) a(2) W a(4)];
      gAxesPosition = [a(1)+os+W a(2) a(3)-os-W a(4)];
   end
end;
return

function h = gcda(hfig, haxes)
%GCDA Get current data axes

h = datachildren(hfig);
if isempty(h) || any(h == haxes)
  h = haxes;
else
  h = h(1);
end
return

function updateLegend(axesHandle)
%UPDATELEGEND Updates the legend if any, in the axes specified by axesHandle

legH = legend(axesHandle);
if ~isempty(legH) && ishandle(legH)
   %legend('RecordSize',axesHandle);
   legend(legH) % Update legend
end
return