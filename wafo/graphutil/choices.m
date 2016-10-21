function choices(name,header,labels,callbacks,inter)
%CHOICES Create a list of choices with uicontrols and callbacks.
%   CHOICES('NAME',HEADER,BUTTONLABELS,CALLBACKS) creates
%   a window with registered name NAME.  The window contains
%   the string HEADER and buttons labeled with BUTTONLABELS.
%   These buttons register callback strings from CALLBACKS.
%   An additional button, labeled 'Close', is added to each
%   choicelist.
%
%   CHOICES is useful for constructing demo menus.
%   Use CHOICES in conjunction with CHOICEX, as in
%   this example.
%       header = 'Easy Example';
%       labels = str2mat('Choice 1','Choice 2','Choice 3');
%       callbacks = str2mat('image(magic(1))','image(magic(2))', ...
%           'image(magic(3))');
%       choices('EXAMPLE', header, labels, callbacks);
%
%   A final, optional, non-zero argument for CHOICES causes the
%   buttons to have their interruptible property set to 'on'.  This
%   is necessary for any demos requiring any mouse button presses, for
%   selections of options or data, for example.

%%%   Loren Shure, 8-14-92.
%%%   Copyright 1984-2000 The MathWorks, Inc.
%%%   $Revision: 5.16 $  $Date: 2000/06/01 03:46:32 $

global CHOICELIST
global CHOICEHANDLES
c = computer;
if ~ischar(name) || ~ischar(header) || ~ischar(labels) || ~ischar(callbacks)
   error('Requires string arguments.');
end
if nargin < 4
   error('Not enough input arguments.')
end
if nargin == 4
   inter = 0;
end
if inter
   yn = 'on';
else
   yn = 'off';
end
uicok = strcmp(c(1:2),'PC') | strcmp(c(1:2),'MA');
if isunix || ~uicok
   uicok = strcmpi(get(0,'TerminalProtocol'),'x');
end
%can't use uicontrols -use menu stuff instead- this is for terminals -UNIX & VMS
if ~uicok
   labels = str2mat(labels,'Done');
   nl = size(labels,1);
   % build up menu string for evaluation
   % fix quotes, if there are any
   ss = deblank(labels(1,:));
   ss = ss(sort([1:length(ss) find(ss=='''')]));
   args = ['''',ss,''''];
   header = header(sort([1:length(header) find(header=='''')]));
   for i = 2:nl
      ss = deblank(labels(i,:));
      ss = ss(sort([1:length(ss) find(ss=='''')]));
      args = [args, ',''', ss,''''];
   end
   k = 1;
   while k > 0 && k < nl
      k = eval(['menu(''',header,''',', args,');']);
      if k == nl || k == 0
         return
      else
         ceval(callbacks(k,:));
      end
   end
   return
end
% can use uicontrols
name = deblank(name);
if isempty(name)
   error('Requires non-blank string argument.')
end
% ensure list doesn't go into figure 1
figs = sort(get(0,'Children'));
openfigs = size(figs);
if ~isempty(figs)
   cf = gcf;
   if cf == 1
      cf = [];
   end
else
   cf = [];
end
fig1 = 1;
if isempty(figs)
   CHOICELIST = [];
   CHOICEHANDLES = [];
   figs = figure('visible','off');
   fig1 = 0;
end
if figs(1) ~= 1
   figs = [figure('visible','off'); figs];
   fig1 = 0;
end
matchn = 0;
for i = 1:size(CHOICELIST,1)
   if strcmp(name,deblank(CHOICELIST(i,:)))
      matchn = i;
      break;
   end
end
if ~matchn
   CHOICEHANDLES = [CHOICEHANDLES(:); 0];
   if isempty(CHOICELIST)
      CHOICELIST = name;
   else
      CHOICELIST = str2mat(CHOICELIST, name);
   end
   matchn = size(CHOICEHANDLES,1);
else
   matchh = 0;
   for i = 1:size(figs,1)
       if figs(i) == CHOICEHANDLES(matchn)
          matchh = i;
          break;
       end
   end
   if matchh
       figure(CHOICEHANDLES(matchn));
       return
   end
end
ss = get(0,'ScreenSize');
xedge = 30;
ybord = 30;
width = 30;
yedge = 35;

% Determine size of text labels

ha = axes('visible','off');
%hh = text(.1,.1,labels,'units','pixel')
hh = text(ones(size(labels,1)+1,1),ones(size(labels,1)+1,1),str2mat(labels,'Close'),'units','pixel');
maxwidth = 0;
height = 0;
for i = 1:length(hh),
    ext = get(hh(i),'extent');
    maxwidth = max(maxwidth,ext(3));
    height = max(height,ext(4));
end      
delete(hh);delete(ha);
yedge = 1.5*height;
height = 6*yedge/7;

imax = 1;

twidth = maxwidth;
% now figure out total dimensions needed so things can get placed in pixels
mwwidth = twidth + width + 2*xedge;
mwheight = (size(labels,1)+2.5)*yedge;
swidth = ss(3); sheight = ss(4);
left = 20;
bottom = sheight-mwheight-ybord;
rect = [left bottom mwwidth mwheight];
CHOICEHANDLES(matchn) = figure('Position',rect,'number','off', ...
       'name','','resize','off','colormap',[],...
       'Menubar','none','visible','off');

fg = CHOICEHANDLES(matchn);
fgs = CHOICEHANDLES(CHOICEHANDLES ~= fg);
set(fgs,'visible','off')
set(gca,'Position',[0 0 1 1]); axis off;

% Place header

hdrpos = [.05 1-1/(size(labels,1)+1.6) .9 1/(size(labels,1)+1.6)];

hh=uicontrol(fg,'style','text','units','normal',...
           'position',hdrpos,'string',header,...
           'Horizontal','center');
set(hh,'backg',get(gcf,'color'),'foreg',[1 1 1]-get(gcf,'color'))
% set up pre-amble so figure 1 is available for rendering in
sb = ['figure(1),set(1,''visible'',''on'');set(',int2str(fg),',''visible'',''off'');'];
se = ';global CHOICEHANDLES;set(CHOICEHANDLES(length(CHOICEHANDLES)),''visible'',''on''),clear CHOICEHANDLES';
for ii=size(labels,1):-1:1
    i = size(labels,1) + 2 - ii;
    h1 = uicontrol(fg,'position',[xedge  (i-.5)*yedge width+twidth height],...
         'callback',[sb callbacks(ii,:) se],...
         'string',['  ',deblank(labels(ii,:)), '  '],...
         'HorizontalAlignment','left','interruptible',yn);
% left justify string inside button
end
% Create Close button
uicontrol(fg,'position',[xedge .5*yedge width+twidth height],...
      'string','  Close  ',...
      'callback',['choicex(''',name,''')'],...
      'HorizontalAlignment','left');

set(fg,'HandleVisibility','callback','visible','on')

if ~isempty(cf)
   figure(cf)
else
   if ~fig1
      close(1);
   end
end

