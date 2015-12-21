function recfig12
% RECFIG12 Probability of exceeding V: Model (dash); data (dots)
% 


% History
% Revised pab 12.06.2001
%  changed from edf(V/Vrms,0,[x F],2); to  plotedf(V/Vrms,[x F],2);
% By pab

% TODO % Fix error in call to cdfmargcnd2d

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end

global V Vrms sphat
x=linspace(0,3.5)';
F=cdfmargcnd2d(x,20,sphat,0);
plotedf(V/Vrms,[x F],2);
hold on,
semilogy(x(F<1),1-F(F(:)<1),'r--')
hold off
xlabel('v')
ylabel('1-F(v)')
title('')
axis([0 3.5 1e-4 1])
axis square
wafostamp('Figure 12','(NR)')