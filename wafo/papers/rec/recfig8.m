function recfig8
% RECFIG8 Estimated Weibull shape parameter versus h=H/Hrms (circles)
%         and 95% pointwise confidence intervals (dash) compared
%         with the final smoothed estimate (solid).
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global phat sphat 

plotmargcnd2dfit(phat,sphat,2)
ylabel('B(h)')
xlabel('h')
wafostamp('Figure 8','(NR)')
