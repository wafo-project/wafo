function wafofig9
% WAFOFIG9  Intensity of rainflow cycles computed from St
%           through Markov approximation, compared with the 
%           cycles found in the simulation.
%

global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

% Use the speed and nit as for the North Sea data
global St fmm xt NNp Nnit Nspeed

paramu = [-6 6 20];
% Only need to calculate Globals which is not empty
% Poseidon / Japan Sea data
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

f.f= mctp2rfm({fmm.f , []});
f.cl=qlevels(f.f);
tp = dat2tp(xt);
rfc = tp2rfc(tp);
if 0,
  cocc(paramu,rfc,f.f)
else
  ccplot(fliplr(rfc));hold on
  pdfplot(f,'k-'),  hold off
end
wafostamp('Figure 9','(NR)')    

return
