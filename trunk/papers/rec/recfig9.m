function recfig9
% RECFIG9 Conditional mean (solid, circle) and standard deviation (dash,cross) of V given H:
%         fitted model (solid, dash); reconstructed field data (circle, cross)       
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global V H Hrms Vrms sphat res 
plotmargcnd2dmom(V/Vrms,H/Hrms,sphat,2,res);
wafostamp('Figure 9','(NR)')
