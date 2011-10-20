function plotmarg2dcdf(V,H,varargin) 
%PLOTMARG2DCDF Plot conditional CDF of X1 given X2=x2
%               and optionally compares it with distribution defined
%               by the parameter vector phat.
%  
%  CALL: plotmarg2dcdf(x1,x2,phat,res,plotflag,figdata,sym);
% 
%       x1,x2 = data
%       phat  =  parameter structure array
%       res   = resolution (default range(x2)/12)
%    plotflag = 0  no plotting
%               1 plot cdf F(x1|x2)                  (default)
%               2 plot 1-F(x1|x2) on a semilog y-scale 
%               3 plot F(x1|x2) on a log log scale
%     figdata = [rows cols Nfig], gives number of rows, columns of which
%               each figure is divided into and and the total number of figures.
%               (default [3 3 1])
%         sym = {s1,s2} cell array of plot symbols for 
%               plotting empirical and theoretical cdf, respectively.
%               (default {'b.','r--'})
%  
%  PLOTMARG2DCDF plots the  empirical CDF of X1 given X2 or 
%  X2 given X1 and compares it a 2D Weibull distribution with parameters
%  given by phat.
% 
% NOTE:  SYM can be given anywhere after X1 and X2
% 
% Example:
%  R = rndray(2,1000,2); x1 = linspace(0,10)';
%  phat = createfdata('distribution',@pdfmarg2d,'params',[2 2 .5])
%  phat.pdfoptions.distribution={'pdfray','pdfray'};
%  phat.pdfoptions.numpar =ones(1,2);
%  d1 = data_1d();
%  survivalPlotflag = 10;
%  plotflag = d1.plotscaleflag('ylog')+ d1.plottypeflag('plot')+survivalPlotflag;
%  plotmarg2dcdf(R(:,1),R(:,2),phat,[],plotflag,[3 3 1],{'k-','g-'});
%
%  See also  cdfmarg2d, edf

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


%  tested on: matlab 5.2
% history
% by pab 25.10.2000


% default values
%~~~~~~~~~~~~~~~
sym = {'b.','r--'}; % default plot symbols for the empirical
                              %  theoretical pdf,respectively
phat =[];
res    = [];
condon = 2;
flag   = 1;
row=3;col=3;Nfig=1;
figdata = [];
ih = ishold; % save hold state



P  = varargin;
Np = length(P);
if Np>0
  strix = zeros(1,Np);
  cellix = strix;
  for ix=1:Np, % finding symbol strings or cell array of symbol strings
    strix(ix)  = ischar(P{ix});
    cellix(ix) = iscell(P{ix});
  end
  k  = find(strix);
  k1 = find(cellix);
  if any(k)
    Nk = length(k);
    if Nk>2,  warning('WAFO:PLOTMARG2DCDF','More than 2 strings are not allowed'),    end
    iy = 1;
    for ix=k      
      sym{iy} = P{ix};
      iy=iy+1;
    end
    Np = Np-Nk;
    P  = {P{find(~strix)}}; % remove strings from input
  elseif any(k1) % cell array of strings
    tmp = P{k1};
    Nk = length(tmp);
    if Nk>2,  warning('WAFO:PLOTMARG2DCDF','More than 2 strings are not allowed'),    end
    iy = 1;
    for ix=1:min(Nk,2)
      if ~isempty(tmp{ix}) && ischar(tmp{ix}), sym{ix}=tmp{ix};end
    end
    Np = Np-1;
    P  = {P{find(~cellix)}}; % remove cell array of strings from input
  end
  if Np>0 && ~isempty(P{1}), phat   = P{1};end
  if Np>1 && ~isempty(P{2}), res    = P{2};end
  if Np>2 && ~isempty(P{3}), flag    = P{3};end
  if Np>3 && ~isempty(P{4}), figdata = P{4};end
end

if isempty(res)
  if condon==1,
    res=range(V(:))/12;
  else
    res=range(H(:))/12;
  end
end

if flag<1,  return,end

nf=length(figdata);
if (nf>0) 
  if ~isnan(figdata(1)),         row  = figdata(1);end
  if (nf>1) && ~isnan(figdata(2)),col  = figdata(2);end
  if (nf>2) && ~isnan(figdata(3)),Nfig = figdata(3);end
end

Nmesh=40;
v1=[];cdfgH=[];
if condon==2,
  
  Xc      = V;
  grp     = floor(H/res)+1; % dividing the data into groups 
  Ngrp    = max(grp);
  h1      = linspace(res/2, (Ngrp-0.5)*res, Ngrp)';
  if ~isempty(phat)
    v1      =linspace(eps,max(V)+range(V)/4,Nmesh);
    [X1,X2] = meshgrid(v1,h1);
    cdfgH   = cdfmarg2d(X1,X2,phat,'condon',condon);
  end
  %max(cdfgH')
  %min(cdfgH')
  xmax    = max(V);
 
else
  Xc      = H;
  grp     = floor(V/res)+1; % dividing the data into groups 
   Ngrp    = max(grp);
  h1      = linspace(res/2, (Ngrp-0.5)*res, Ngrp)';
  if ~isempty(phat)
    v1      = linspace(eps,max(H)+range(H)/4,Nmesh);
    [X1,X2] = meshgrid(v1,h1);
    cdfgH   = cdfmarg2d(X2,X1,phat,'condon',condon);
  end
  xmax    = max(H);
end


%xmax=min(xmax,4);

fignr    = gcf;
fignrold = fignr;
%figure(fignr)

iy=0;
for ix=Ngrp:-1:min(grp),
  if iy==row*col,iy=0; 
    fignr=fignr+1;
    if Nfig<=fignr-fignrold, break, end
    figure(fignr);
  end
  tmp = Xc(grp==ix);%find data in group number ix
  if length(tmp)>max(3,0),% if less than 6 observations in the group 
    iy=iy+1;
    subplot(row,col,iy)
    if ih, hold on, end % make sure we have hold on for each subplot
    
    Femp = edfcnd(tmp,0,[],'wdata',true);
    Ns = 2;
    Femp.labels = {'x1',['F(x1|x2=' num2str(h1(ix),Ns) ')']};
    plot(Femp,flag,sym{1})
    if ~isempty(phat)
       Fth = wdata(cdfgH(ix,:)',v1');
       Fth.labels = Femp.labels;
       hold('on'), 
       plot(Fth,flag,sym{2})
       hold('off')
    end
    %grid on
    
    axis square

  end
end

if ~ih, hold off, end
