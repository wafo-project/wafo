function plothmm(x,z,t,z_range,Title,Ylabel,colType,fontsize)
% PLOTHMM  plots a Hidden Markov Model.
%
% CALL: plothmm(x,z)
%       plothmm(x,z,t,z_range,Title,Ylabel,colType,fontsize)
%
%   x       = Switching process.
%   z       = Regime process.
%   t       = Time.
%   z_range = Range of regime, e.g. [1 3] for 3 regime states.
%   Title   = Title of plot.
%   Ylabel  = Y-label of x-plot.
%   colType = 0: One colour (default), 1: Different colours for each regime
%   fontsize= Fontsize of text.
%
% Examples: Switching AR(1)-process. (Example 3 in thesis)
%   P= [0.975 0.02 0.005; 0.01 0.98 0.01; 0.005 0.02 0.975];
%   C = [1 1 1]';
%   A = [1 -0.5; 1 -0.3; 1 0.5];
%   m = [-1 0 3]';
%   s2 = [1 1 1.44]';
%   [x,z] = sarmasim(C,A,m,s2,P,500);
%   plothmm(x,z)
%   plothmm(x,z,(0:499)/400,[1 3],'Switching AR(1)-process','X(t)')
%   plothmm(x,z,[],[1 3],'','',1,20)

% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

if nargin<3, t=[]; end
if nargin<4, z_range=[]; end
if nargin<5, Title = ''; end
if nargin<6, Ylabel=''; end
if nargin<7, colType=0; end
if nargin<8, fontsize=[]; end

if isempty(t)
  t = (1:length(x))';
end

if isempty(z_range)
  zmin = min(z); zmax = max(z);
else
  zmin = z_range(1); zmax = z_range(2);
end

tmin = min(t);
tmax = max(t);

subplot(2,1,1)
if colType == 0
  plot(t,x)
elseif colType == 1
  col = get(gca,'ColorOrder');
  i0=1; n=length(x);
  slut=0;
  while ~slut
    I=(find(z(i0:n)~=z(i0)));
    if ~isempty(I)
      i1 = min(I) + i0-1;
    else
      i1=n; slut = 1;
    end
    plot(t(i0:i1),x(i0:i1),'Color',col(z(i0),:));
    i0=i1;
    hold on
  end
end

v=axis; axis([tmin tmax v(3:4)])
h1=gca;
if ~isempty(fontsize), set(h1,'Fontsize',fontsize); end
pos1 = get(h1,'Position');
title(Title); ylabel(Ylabel);

subplot(2,1,2)
if colType == 0
  stairs(t,z)
elseif colType == 1
  i0=1; n=length(x);
  slut=0;
  while ~slut
    I=(find(z(i0:n)~=z(i0)));
    if ~isempty(I)
      i1 = min(I) + i0-1;
    else
      i1=n; slut = 1;
    end
    [XX,YY] = stairs(t(i0:i1),z(i0:i1));
    plot(XX,YY,'Color',col(z(i0),:));
    i0=i1;
    hold on
  end
end

axis([tmin tmax zmin-0.5 zmax+0.5]);
h2=gca;
if ~isempty(fontsize), set(h2,'Fontsize',fontsize); end

pos2 = get(h2,'Position');

Pos2 = pos2; Pos2(4)=.1;
Pos1 = pos1; Pos1(2) = Pos2(2)+Pos2(4)+0.02; Pos1(4) = pos1(4)+pos1(2)-Pos1(2);

set(h1,'Position',Pos1);
set(h2,'Position',Pos2);

set(h1,'XTickLabel',[]);
set(h2,'YTick',zmin:zmax);
