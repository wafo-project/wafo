function recfig10
% RECFIG10 Joint distribution of  V and H:
%         Fitted distribution (solid); Kernel density (dash); data (dots)
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global V H Hrms Vrms ft2 fkde 

plot(V/Vrms,H/Hrms,'k.')
axis([0 3.5 0 3.5]), axis square
hold on
pdfplot(ft2,1,'k')
%hold on
pdfplot(fkde,1,'k--')
title('Joint distribution of v and h')
hold off
wafostamp('Figure 10','(NR)')
