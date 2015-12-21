function sphat=margcnd2dsmfun2(phat,h2,csm,lin,visual,f0)
%MARGCND2DSMFUN2 Smooths the MARGCND2D distribution parameters. 
% 
% CALL:  sphat = dist2dgsmfun2(phat,x2,csm,lin)
% 
%    sphat = smoothed parameter structure array evaluated at x2
%    phat  = parameter structure array
%    x2    = evaluation points (default phat.x{1}(:,1)) 
%    csm   = [csma csmb csmc], smoothing parameter vector which defines
%            the smoothing of parameter A,B and C, respectively 
%            0 -> LS straight line
%            1 -> cubic spline interpolant (default [1 1 1])
%    lin   = [lina linb linc], vector defining the extrapolation of
%            the parameters A,B and C,  respectively (default [1 1 1])
%            0 : No linear extrapolation outside the range of data
%            1 : Linear extrapolation outside the range of data 
%  visual  = 0 : No visual fitting by the user (default)
%            1 : Visual fitting by the user by specifying points using
%                the mouse and interpolate the points with a spline afterwards.
%
%  See also   fitmargcnd2d cssmooth 

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Adapted to  cssmooth  by GL Feb 2011
% Tested on: Matlab 5.2
% history
% revised pab July2004  
% revised pab 20.01.2004  
% revised pab 08.02.2000
%  - added visual
% by  Per A. Brodtkorb 28.10.98


if nargin<2 ||isempty(h2),
  ind0 = (~isnan(phat.x{1}(:,1)));
  h2=phat.x{1}(ind0,1);
end
h2=h2(:);

if nargin<3||isempty(csm)
  csm= [];
end
if (nargin<4) || isempty(lin),
    lin=[];
end
if (nargin<5) || isempty(visual),
    visual=0;
end
if (nargin<6) 
  f0 = [];
end
sphat=phat;

if isfield(phat,'CI')
  sphat=rmfield(sphat,'CI');
end
sphat.date=datestr(now);
pvhsiz=size(phat.x{1});
 
n=length(h2(:));
spvhsiz=[n pvhsiz(2)];
sphat.x{1}=zeros(spvhsiz);
sphat.x{1}(:,1)=h2;

useSTATS = (1 & isfield(phat,'stats1'));
if useSTATS && ~isempty(phat.stats1) && strncmpi(phat.dist{1},'gamma',2)
   CSMA = 1;
   CSMB = 1;
   linA = 1;
   linB = 1;
   if (~isempty(csm) && ~isempty(csm)), CSMA = csm(1); end
   if (~isempty(csm) && length(csm)>1), CSMB = csm(2); end
   if (~isempty(lin) && ~isempty(lin)), linA = lin(1); end
   if (~isempty(lin) && length(lin)>1), linB = lin(2); end
    %smooth on conditional mean and variance
    ind1   = find(~isnan(sum([phat.stats1{:}],2)));
    hi = phat.stats1{1}(ind1);
    
    
    mi = phat.stats1{2}(ind1);
    vi = phat.stats1{3}(ind1);
    Nptsi = phat.stats1{4}(ind1);
    dab =  hi.^2./(Nptsi);
    %h2 = unique(h2);
    Mv = cssmooth(hi,mi,CSMA,h2,linA,dab);
    if any(Mv<=0) 
      dab1 =  log(mi)+log(mi+dab);
       Mv = exp(cssmooth(hi,log(mi),CSMA,h2,linA,dab1));
    end
    Sv = cssmooth(hi,sqrt(vi),CSMB,h2,linB,dab);
    if any(Sv<=0) 
      dab1 =  log(sqrt(vi))+log(sqrt(vi)+dab);
       Sv = exp(cssmooth(hi,log(sqrt(vi)),CSMB,h2,linB,dab1));
    end
    %plot(hi,mi,'.',hi,sqrt(vi),'.',h2,Mv,h2,Sv)
    %drawnow
    %pause
    if 1
      Av = Mv.^2./Sv.^2;
      Bv = Sv.^2./Mv;
    else
      Av = cssmooth(h2,Mv.^2./Sv.^2,0.99,h2);
      Bv = cssmooth(h2,Sv.^2./Mv,0.99,h2);
    end
    Cv = 0;
    sphat.x{1}(:,2) = Av;
else
  [sphat.x{1}(:,2), Bv , Cv]=margcnd2dsmfun(phat,h2,csm,lin);
end
  if ~isempty(Bv)
    sphat.x{1}(:,3)=Bv;  
  end

  if ~isempty(Cv) && any(Cv~=0),
    sphat.x{1}(:,4)=Cv;
    
  end
sphat.csm=csm;
sphat.lin=lin;

if ~visual, return, end

% Visual Fitting of the parameters starts here:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% this may be as good as any other method when extrapolating the
% parameter outside the observed range.

sphat0 = sphat;

fig0 = gcf;
figure(fig0)
x1 = linspace(0,5);
clf
cmpDistributions(x1,sphat,f0)
figure(fig0+1)
plotmargcnd2dfit(phat,sphat)


for ix = 1:10
  figure(fig0+2)
  plotmargcnd2dfit(phat,sphat,1)
  if (ix>1),
    hold on
    plotmargcnd2dfit(phat,sphat0,1)
    hold off
  end
  [x,y]=graphicinput;
  
  sphat.visual=~isempty(x);
  if ~isempty(x),
    sphat.x{1}(:,2)=cssmooth(x,y,1,sphat.x{1}(:,1),1);
    hold on, plot(h2,sphat.x{1}(:,2),'g-'),hold off % plot visual fit
    figure(fig0)
    cmpDistributions(x1,sphat,f0);
    figure(fig0+1)
    plotmargcnd2dfit(phat,sphat)
  end

  if ~isempty(Bv)
    figure(fig0+3)
    plotmargcnd2dfit(phat,sphat,2)
    if (ix>1),
      hold on
      plotmargcnd2dfit(phat,sphat0,2)
      hold off
    end
    [x,y]=graphicinput;
    sphat.visual(2)=~isempty(x);
    if ~isempty(x),
      sphat.x{1}(:,3)=cssmooth(x,y,1,sphat.x{1}(:,1),1);
      hold on, plot(h2,sphat.x{1}(:,3),'g-'),hold off % plot visual fit
      figure(fig0)
      cmpDistributions(x1,sphat,f0)
      figure(fig0+1)
      plotmargcnd2dfit(phat,sphat)
    end
  end
  if ~isempty(Cv) && any(Cv~=0),
    if strcmpi(phat.dist{1}(1:2),'ra'), col=2;else col=3;end
    figure(fig0+4)
    plotmargcnd2dfit(phat,sphat,col)
    if (ix>1),
      hold on
      plotmargcnd2dfit(phat,sphat0,col)
      hold off
    end
    [x,y]=graphicinput;
    sphat.visual(col)=~isempty(x);
    if ~isempty(x),
      sphat.x{1}(:,col+1)=cssmooth(x,y,1,sphat.x{1}(:,1),1);
      hold on, plot(h2,sphat.x{1}(:,2),'g-'),hold off % plot visual fit
      figure(fig0)
      cmpDistributions(x1,sphat,f0)
      figure(fig0+1)
      plotmargcnd2dfit(phat,sphat)
    end 
  end
  if ix<10
    ButtonName=questdlg('Do you want to visually fit once more?', ...
                       'Yes','No');
    if ~strncmpi(ButtonName,'yes',3)
      break
    end
  end
end

return


function [x, y]= graphicinput
% Local function which reads points selected with the mouse from screen
% 
  x=[];y=[];
  n=0;
  if 0
  lim=axis; % find the current scaling of the figure
  %answer = input(['Enter scaling for the current plot. (default [' num2str(lim,4) '] )'])
  prompt={'Enter scaling for the current plot.'};
  answer=inputdlg(prompt,'Visual fitting',1,{ num2str(lim,4) },'on');	
  Na=length(answer);
  if Na>0,
    answer{1}=str2num(answer{1});
    if (isempty(answer{1})||all(answer{1}==lim)),  else axis(answer{1}), end
  end
  end
  disp('Left mouse button picks points')
  disp('Middle mouse button (or z) to zoom in/out')
  disp('Right mouse button picks last point')
  disp('(if less than 6 points are selected then no fitting is made.)')
  
  leftButton    = 1;
  middleButton  = 2;
  button = 1;

  hold on
  while (~isempty(button) && (button==leftButton)),
    [xi,yi,button] = ginput(1);
    if ~isempty(button)
      if (button == leftButton),
	plot(xi,yi,'go')
	n=n+1;
	x(n,1)=xi;
	y(n,1)=yi;
      elseif( (button == middleButton) || strcmpi(char(button),'z') )
	disp('click the left mouse button to zoom in on the point')
	disp('under the mouse.')
	disp('Click the right mouse button to zoom out.')
	disp('Click and drag to zoom into an area.')
	%disp('To resume the identification of points:)
	zoom on
	button = input('Hit the ENTER key on keyboard when finished zooming.');
	zoom off
	button = leftButton;
      end
    end
  end
  
  if n<6,
    x=[];y=[];
  end
  hold off

return

function cmpDistributions(x1,sphat,f0)
 f = pdfmargcnd2d(x1,x1,sphat,'wdata',true','meshgrid',true);
 plot(f);
 if ~isempty(f0), 
   hold on
   pdfplot(f0,'r')
   hold off
 end
 return
