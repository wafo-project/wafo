function wafofig2
% WAFOFIG2  Spectral densities:
%           Original Torsethaugen spectrum (solid),
%           Estimated from linear simulations (dash), and
%           95% Confidence interval (dotted)
%

global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

global St Ste
plotspec(St,1,'k-')
ih=ishold;
hold on
plotspec(Ste,1,'--')
if ~ih, hold off, end
wafostamp('Figure 2','(ER)')    
