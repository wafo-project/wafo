%KDEDEMO1 Demonstrate the smoothing parameter impact on KDE
%
% KDEDEMO1 shows the true density (dotted) compared to KDE based on 7
% observations (solid) and their individual kernels (dashed) for 3
% different values of the smoothing parameter, hs.
%
% Example
% kdedemo1;
%
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -*- Mode: Matlab -*- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kdedemo1.m --- 
%% Author          : Per Andreas Brodtkorb
%% Created On      : Sat Mar 06 10:54:08 2004
%% Last Modified By: Per Andreas Brodtkorb
%% Last Modified On: Sat Feb 05 10:53:14 2005
%% Update Count    : 66
%% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 x      = linspace(-4,4);
 x0     = x/2;
 data   = rndnorm(0,1,7,1);
 kernel = 'gaus';
 hs     = hns(data);
 % hs     = hldpi(data,kernel);
 hVec   = [hs/2, hs, 2*hs];
   
 kopt = kdeoptset('kernel','gaus');
 for ix=1:3 
   figure(ix)
   h = hVec(ix);
   kopt.hs = h;
   f2 = kde(data,kopt,x);
   pdfplot(f2,'k-');
   title(sprintf('h_s = %s', num2str(h,2)));
   ylabel('Density');
   
   hold on;
   plot(x,pdfnorm(x,0,1),'k:');
   n = length(data);
   plot(data,zeros(size(data)),'bx','Markersize',10);
   y = mkernel(x0,kernel)/(n*h); 
   for i=1:n
     plot(data(i)+x0*h,y,'b--');
     plot([data(i) data(i)], [0 max(y)],'b');
   end
   set(gca,'ytick',0:.1:0.5);
   axis([min(x),max(x), 0 0.5]);
   if ismatlab,
    axis fill;
   else,
    axis auto;
   end
   axis([min(x),max(x), 0 0.5]);
   hold off;
   if (ix==3)
     xlabel('x');
   end
   %exportfig(gcf,sprintf('kdedemo1f%d.eps',ix),'height',2.5,'width',5);
 end
