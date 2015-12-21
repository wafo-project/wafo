function black(lty)
%BLACK  Presets for making b/w graphs
%
% CALL:   black(lty)
%         black off
%
% lty    = optional string of linestyles separated by |.
%          (default '-|--|-.|:|:.|:o|:x|:+|:*|:s')
%
% BLACK sets ColorOrder and LinestyleOrder of current figure with run
% through of the different linestyles instead of colors.   
% (More information found at
% Helpdesk->HandleGraphicsObjectProperties->Axes->LineStyleOrder.)
%
%Example
% t=linspace(0,4*2*pi)'; x = sin(t); y = cos(t);
% figure(1),plot(t,x,t,y), figure(2), black,plot(t,x,t,y)
% black off % reset ColorOrder and LinstyleOrder to factory settings 
%
% See also: plot 


% revised pab 
% - changed help header to WAFO style
% - added off
%Time-stamp:<Last updated on 00/07/12 at 09:31:30 by even@gfi.uib.no>
%File:<d:/home/matlab/black.m>

error(nargchk(0,1,nargin));
if nargin < 1 || isempty(lty) || ~ischar(lty)
  lty='-|--|-.|:|:.|:o|:x|:+|:*|:s';   %default LineStyleOrder
end 
CO = [0 0 0]; %Default ColorOrder is black
if strcmpi(lty,'off'),
  %ltyold = get(gcf,'DefaultAxesLineStyleOrder');
  %COold  = get(gcf,'DefaultAxesColorOrder');
  lty = '-';
  CO  = [0      0     1;...
        0      0.5   0;...
        1      0     0;...
        0      0.75  0.75;...
        0.75   0     0.75;...
        0.75   0.75  0;...
        0.25   0.25  0.25;];
end
set(gcf,'DefaultAxesLineStyleOrder',lty);
set(gcf,'DefaultAxesColorOrder',CO);

% Avoid setting nextplot to add, since this makes MATLAB remember all
% previous lines regardless of new plots in figure window (???)
%set(gca,'nextplot','add','colororder',[0 0 0],'linestyleorder','-|--|:|-.')






