function recfig6
% RECFIG6 Variability of simulated e(g(u)-u) (circles)
%         for Gaussian processes with the same spectral densities and
%         lengths as the reconstructed data set compared with the
%         observed e(g(u)-u) for the reconstructed data set (dash)
% 

% revised pab 2 Sept 2005
%  -updated order of arguments to mctrtest once again.
% revised pab 18.04.2001
% updated call to mctrtest.

global RECFIGNUM
if isempty(RECFIGNUM)
  disp('You must start recdemo in order to run this script')
  return
end

global Sr test0 csm1 csm2 param 

%Ns=100; %This takes a while and is used in the article
Ns=50;
opt = troptset('csm',csm1,'gsm',csm2,'param',param);
test1 = testgaussian(Sr,[36000,Ns],test0,[],opt);
wafostamp('Figure 6','(NR)')
