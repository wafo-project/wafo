function recfig11
% RECFIG11 Probability of exceeding H: Model (dash); data (dots)
% 

% Tested on: Matlab 5.x
% History:
% Revised pab 12.06.2001
%  changed from edf(H/Hrms,0,[x F],2); to  plotedf(H/Hrms,[x F],2);

% TODO % Fix error in call to cdfmargcnd2d
global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global H Hrms sphat 
x=linspace(0,3.5)';
F=cdfmargcnd2d(20,x,sphat,0);
plotedf(H/Hrms,[x F],2);
xlabel('h')
ylabel('1-F(h)')
title('')
hold on,
semilogy(x(F<1),1-F(F(:)<1),'r--')
hold off
axis([0 3.5 1e-4 1])
axis square
wafostamp('Figure 11','(NR)')
