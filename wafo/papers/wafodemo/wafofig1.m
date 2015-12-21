function wafofig1
% WAFOFIG1  Linear simulation from a Torsethaugen spectrum.
%           Use the mouse to zoom in on a specific region.
%

global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

global xt Hm0 Tp

disp(['        Hm0 = ' num2str(Hm0), '   Tp = ', num2str(Tp)])
waveplot(xt,'-',1,1)
zoom xon
wafostamp('Figure 1','(ER)')    
