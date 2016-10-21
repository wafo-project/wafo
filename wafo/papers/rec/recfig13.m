function recfig13
% RECFIG13 The conditional probability of exceeding V given H: Model (dash); data (dots)
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global H Hrms V Vrms sphat res 



plotmargcnd2dcdf(V/Vrms,H/Hrms,sphat,2,res,2);
wafostamp('Figure 13','(NR)')