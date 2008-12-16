function wafofig10
% WAFOFIG10 Intensity of trough-crest cycles computed from  St
%           through Markov approximation, compared with the 
%           cycles found in the simulation.
%

global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end
% TODO % Is not finished
global St fmm xt NNp Nnit Nspeed

paramu = [-6 6 20];
% 
if isempty(fmm)
  disp('This takes several  minutes to finish => several hours ... ')
  disp('depending on input arguments and your computer')
  % Calculate the theoretical distribution
  if Nnit<0
    opt = rindoptset('method',abs(Nnit),'speed',Nspeed);
  else
    opt = rindoptset('method',0,'nit',(Nnit),'speed',Nspeed);
  end
  fmm = spec2mmtpdf(St,0,'mm',[0 7 NNp], paramu,opt);
end
f=fmm;
if 1,
  % WAT CALL
  f.f = mctp2tc(fmm.f,0,paramu);
  %else
  %% WAFO new call not implemented
  %f.f = mctp2rfm({fmm.f , []});
  %f.f = rfm2tcpdf(f.f,0,paramu);
end
f.cl=qlevels(f.f);
tc  = dat2tc(xt);
[mM Mm] = tp2mm(tc);
if 0,
  cocc(paramu,Mm,f.f)
else
  ccplot(Mm); hold on
  pdfplot(f,'k-'), hold off
end
wafostamp('Figure 10','(NR)')    
return
