function wafofig3
% WAFOFIG3  Probability density distributions (pdf) of wave period, Tt,
%           i.e., down-to-upcrossing:
%           Theoretical pdf given the theoretical spectral density, S (solid), 
%           theoretical pdf given the estimated spectral density from
%           the simulated data, Sest (green dash-dot), and
%           a kernel density estimate from the simulated data (red dash).
%           The histogram shows the wave periods extracted from simulated 
%           data.
%

% Revised pab Feb2005
% -updated call to kdebin and spec2thpdf
global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

global xt u rate
global St Ste Tt fTt fTte Np nit speed  
global kdeTt kernel hs L2

if nit<0
  opt = rindoptset('method',abs(nit),'speed',speed);
else
  opt = rindoptset('method',0,'nit',(nit),'speed',speed);
end

% Only need to calculate Globals which is not empty
if isempty(Tt)
  Tt=dat2wa(xt,u,'d2u','dw',rate);
end
if isempty(fTt)
  fTt = spec2thpdf(St,u,'Tt',[0 11 Np],[],opt);
end
if isempty(fTte)
  fTte = spec2thpdf(Ste,u,'Tt',[0 11 Np],[],opt);
end
if isempty(kdeTt)
  kopt = kdeoptset('kernel',kernel,'hs',hs,'L2',L2);
  kdeTt=kdebin(Tt,kopt);
end
histgrm(Tt,max(22,2*sqrt(length(Tt))),[],1)
hold on
pdfplot(kdeTt,'r--')
pdfplot(fTt,'-')
pdfplot(fTte,'g-.')
hold off
axis([0 inf 0 inf])

wafostamp('Figure 3','(ER)')    