function wafofig8
%WAFOFIG8  Joint distribution (pdf) of crest wavelength, Lc, and crest amplitude, Ac for extremal waves 
%           of a JONSWAP spectrum which is truncated outside 0.5*wp and 1.5wp 
%

% TODO % This is not implemented. Finish the code!
global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

disp('Not implemented yet')
return


% Use the speed and nit as for the North Sea data
global fLcAc NNp Nh Nnit Nspeed 

if Nnit<0
  opt = rindoptset('method',abs(Nnit),'speed',Nspeed);
else
  opt = rindoptset('method',0,'nit',(Nnit),'speed',Nspeed);
end

% Only need to calculate Globals which is not empty
% Poseidon / Japan Sea data
if isempty(fLcAc)
  disp('This takes several  minutes to finish => several hours ... ')
  disp('depending on input arguments and your computer')
  Sj=jonswap;
  % Find the peak frequency
  ind = wfindpeaks(Sj.S);
  % Truncate the spectrum outside 0.5wp and 1.5wp
  Sj.S(1:floor(ind(1)*.5))=0;
  Sj.S(floor(ind(1)*1.5):end)=0;
  Sk=spec2spec(Sj,'k1d');
  % Calculate the theoretical distribution
  fLcAc = spec2thpdf(Sk,0,'LcAc',[0 250 NNp],Nh,opt);
end
pdfplot(fLcAc,'k-')
wafostamp('Figure 7','(NR)')    
