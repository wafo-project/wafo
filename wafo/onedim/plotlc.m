function h = plotlc(lc,plotflag,ma,sa)
%PLOTLC Plots level-crossing spectrum (lc) 
%
% CALL:  h = plotlc(lc,plotflag,ma,sa);
%
%        h = handle of the graphics object 
%       lc = two column matrix with levels and number of upcrossings,
%            i.e., level-crossing spectrum (LCS).
% plotflag = 1  plots LCS
%            2  plots LCS and the theoretical one for a Gaussian
%               process (default) 
%            11 plot level-crossing intensity
%            12 plot level-crossing intensity and the theoretical one for
%               a Gaussian process  
%  
%    ma,sa = mean and standard deviation of the process 
%            (default ma = 0, sa estimated from lc)
%
% Example: 
%  x = load('sea.dat');
%  lc = dat2lc(x,0.2); 
%  plotlc(lc) 
%
% See also  dat2lc, mm2lc, pdfnorm, polyfit

% Tested on: Matlab 6.0
% History:
% revised pab Feb2004  
%  - added plotflag ==11, 12
% revised by jr 06.07.00. 
% - Division by 2 from 05.07.00 removed (line 42). 
% - sa -> sa^2 (line 66)
% revised by jr 05.07.00. Line 39: division by 2
% revised by pab 11.08.99
% - removed the plot procedure from the old mm2cross function  
% - added the option of overplotting of the theoretical one 
%    for a Gaussian process

error(nargchk(1,4,nargin))
if nargin<2||isempty(plotflag)
 plotflag=2;
end
[cmax icmax]=max(lc(:,2));% maximum crossing
if (nargin <4||isempty(sa))&& mod(plotflag,10)==2, 
  % if stdev of x unknown then estimate it
  cros   = lc;
  indz   = (cros(:,2)==0);
  %lncros = cros; 
  logcros(indz,2)  = inf; % this is done to avoid a warning message
  logcros(~indz,2) = -log(cros(~indz,2));
  logcmin      = logcros(icmax,2);
  logcros(:,2) = sqrt(2*abs(logcros(:,2)-logcmin)); %sqrt(2*logcros)
  logcros(1:icmax,2) = 2*logcros(icmax,2)-logcros(1:icmax,2); 
  p = polyfit(lc(10:end-9,1),logcros(10:end-9,2),1); %least square fit
  
  sa = 1/p(1);% estimated standard deviation of x
  
  %plot(lc(:,1),logcros(:,2)),hold on
  %plot(lc(:,1),polyval(p,lc(:,1)),'r'),hold off, pause
end

if (nargin<3||isempty(ma))&& mod(plotflag,10)==2, 
  ma=0; % the mean is apriori zero.
  %ma=p(2); %, lc(icmax,1) % new estimated mean
  %ma=lc(icmax,1); % maximum cr intensity should be at the mean
  % but this not always the case
end

if plotflag>2
  lc(:,2) = lc(:,2)/cmax;
  cmax = 1;
  ylabtxt = 'Relative crossing intensity';
else
  ylabtxt = 'Number of upcrossings';
end

%  clf
hh = stairs(lc(:,1),lc(:,2));
axis([-1.05*max(abs(lc(:,1))) 1.05*max(abs(lc(:,1))) 0 1.05*cmax])
title('Crossing spectrum')
xlabel('level u')
ylabel(ylabtxt)

% Make a special graph box.
set(gca,'Box','off','TickLength',0.01*[1 1],'TickDir','out')


% bar(d(:,1),d(:,2))
if (mod(plotflag,10)==2)
  y = pdfnorm(lc(:,1),ma,sa.^2);  
  y = cmax*y/y(icmax);%Normalization needed to overplot the crossingspectrum.
  np = get(gca,'NextPlot');    
  set(gca,'NextPlot','add')    
  hh1 = plot(lc(:,1),y,'r--','LineWidth',1);% Plots density line over histogram.
  set(gca,'NextPlot',np) 
end

if nargout == 1 && plotflag==2
  h = [hh; hh1];
elseif nargout==1
  h=hh;
end  
 
return
  
