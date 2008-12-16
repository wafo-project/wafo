function recfig5
% RECFIG5 Transfer function, g, versus the crossing level u
%     of the reconstructed time series estimated by smooting (solid)
%     compared with empirical estimate (irregular) and the theoretical 
%     one for a Gaussian process (dash).
% 

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end
global xr m0 g g2 
trplot(g, g2,mean(xr(:,2)),sqrt(m0))
wafostamp('Figure 5','(NR)')
