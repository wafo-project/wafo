function recfig4
% RECFIG4 Estimated spectral density (solid) and 95% confidence intervals (dots)
%
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end

global Sr
Sf=ttspec(Sr,'f');
plotspec(Sf) 
wafostamp('Figure 4','(NR)')
