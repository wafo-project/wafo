function recfig2
% RECFIG2 10 minutes mean values of wind  (dash) and direction (solid)
%         110 m above mean water level 24-Dec-89 at Statfjord A.
%         Period considered for field data (horizontal solid).  
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end

global wind 
if isempty(wind)
  wind=load('sfa89.dat');
end
  
ih=ishold; 
sym='-';
h=plot([17 20],34*[1 1],sym );
set(h,'linewidth',2.5)
axis([0 24 0 40])
hold on
h=plot([20.333333 21.33333333],34*[1 1],sym );
set(h,'linewidth',2.5)
time= wind(:,1);

switch 1
  case 1, %plotyy
    
    [ax1 h11 h22]=plotyy(time,wind(:,2),time,wind(:,3));
    set(h11, 'LineStyle' , '--')
    xlabel('Time (hours)')
    ylabel('Mean wind speed (m/s)')     
    set(gcf,'currentaxes',ax1(2))
    ylabel('Wind direction (degrees)')
    set(gcf,'currentaxes',ax1(1))
    
  case 2, % quiver plot
    % alternative which maybe is better
    lstep=6; % wind direction every  hour
    theta = -(90+wind(1:lstep:end,3))*pi/180; r = 2*ones(size(theta));
    [u,v] = pol2cart(theta,r);
    h1=plot(time,wind(:,2));hold on
    set(h1,'linewidth',.3) % make thinner line
    % NB: multiply u with 35/24 to avoid axis equal
    h2=quiver(time(1:lstep:end),wind(1:lstep:end,2),u*35/24,v,0,'filled'); %  plot arrows
    set(h2,'linewidth',1) %, set(h2(2),'linewidth',1) % set thicker line
    h2=plot([17 20],ones(1,2)*34,'k-');
    set(h2,'linewidth',1.5)
    h2=plot([20.33333, 21.33333],ones(1,2)*34,'k-'); % period considered
    set(h2,'linewidth',1.5)
    axis([0 25 0 35]),axis square
    xlabel('Time (hours)')
    ylabel('Mean wind speed (m/s)') 
end
grid on

if ~ih,hold off,end

grid on
title('Wind conditions at Statfjord A   24 Dec. 1989')
wafostamp('Figure 2','(ER)')
