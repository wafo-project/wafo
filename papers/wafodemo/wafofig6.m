function wafofig6
% WAFOFIG6  Joint distribution (pdf) of crest front period, Tcf, and crest amplitude, Ac 
%           given that the crest period, Tc=5s:
%           Theoretical joint density (solid), 
%           kernel density estimate (dash) of Tcf and Ac
%           given that 4.8s < Tc < 5.2s of the data from Poseidon in the Japan Sea (dots)
%
%

% Revised pab Feb2005, updated call to kdebin  
  
global  WAFOFIGNUM

if isempty(WAFOFIGNUM)
  disp('You must start wafodemo in order to run this script')
  clear global WAFOFIGNUM
  return
end

global JTcf JAc Jind Jxn Jrate 
global fTcfAcTc JNp Jh Jnit Jspeed 
global kdeTcfAcTc Jkernel Jhs JL2

if Jnit<0
  opt = rindoptset('method',abs(Jnit),'speed',Jspeed);
else
  opt = rindoptset('method',0,'nit',(Jnit),'speed',Jspeed);
end

% Only need to calculate Globals which is not empty
% Poseidon / Japan Sea data

if isempty(JTcf)
  [JVcf, JHd, JAc,JAt,JTcf,JTcr] = dat2steep(Jxn,Jrate,0);
  Jind=find((4.8<JTcf+JTcr).*(JTcf+JTcr<5.2));
end
if isempty(fTcfAcTc)
  disp('This takes several  minutes to finish => several hours ... ')
  disp('depending on input arguments and your computer')
  Sj=dat2spec(Jxn);
  fTcfAcTc = spec2thpdf(Sj,0,'TcfAc',[5 5 JNp],Jh,opt);
end

if isempty( kdeTcfAcTc)
  kopt = kdeoptset('kernel',Jkernel,'hs',Jhs,'L2',JL2);
  kdeTcfAcTc=kdebin([JTcf(Jind) JAc(Jind)],kopt);
  if 1,
    r = evalpdf(kdeTcfAcTc,JTcf(Jind), JAc(Jind),'linear');
    kdeTcfAcTc.cl = qlevels2(r,kdeTcfAcTc.pl); % calculate the levels which encloses fkde.pl
                                               % percent of the data (v,h)
  end
end  
plot( JTcf(Jind), JAc(Jind),'.'), hold on
pdfplot( kdeTcfAcTc,'r--')
pdfplot(fTcfAcTc,'k-')
hold off
wafostamp('Figure 6','(NR)')    