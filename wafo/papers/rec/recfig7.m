function recfig7
% RECFIG7 Estimated Weibull scale parameter versus h=H/Hrms (circles)
%         and 95% pointwise confidence intervals (dash) compared
%         with the final smoothed estimate (solid).
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global phat sphat 

plotmargcnd2dfit(phat,sphat,1)
ylabel('A(h)')
xlabel('h')
wafostamp('Figure 7','(NR)')
