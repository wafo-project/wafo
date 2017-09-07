function [f,Hrms,Vrms] = jhvnlpdf2(Hd,Vcf,Hm0,Tp,gam,normalizedInput)
%JHVNLPDF2 Joint (Vcf,Hd) PDF for non-linear waves with a JONSWAP spectrum.
%
%  CALL: f = jhvnlpdf2(Hd,Vcf,Hm0,Tp,gamma)
% 
%  f     = pdf structure evaluated at meshgrid(Vcf,Hd)
%  Hd    = zero down crossing wave height [m]
%  Vcf   = crest front velocity    [m/s]
%  Hm0   = significant wave height [m]
%  Tp    = Spectral peak period    [s]
%  Gamma = Peakedness parameter of the JONSWAP spectrum
%
% JHVNLPDF2 approximates the joint distribution of (Vcf, Hd), i.e., crest
% front velocity (Ac/Tcf) and wave height, for 2nd order Stokes waves with a
% JONSWAP spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 13 classes of
% GAMMA between 1 and 7. Between 140000 and 150000 zero-downcrossing waves
% were simulated for each class.
% JHVNLPDF2 is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%
% Example:
% Hm0 = 6;Tp = 9; gam=3.5
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,4*2*Hm0/Tp)';
% fl = jhvpdf2(h,v,Hm0,Tp,gam);
% f = jhvnlpdf2(h,v,Hm0,Tp,gam);
% w = linspace(0,40,5*1024+1).';
% S = jonswap(w,[Hm0, Tp, gam]);
% dt = .3;
% x = spec2nlsdat(S,80000,dt); rate = 4;
% [vi,hi] = dat2steep(x,rate,1);
% fk = kdebin([vi,hi],{'L2',.5,'inc',128}); 
% fk.title = f.title; fk.labx = f.labx;  
% plot(vi,hi,'.'), hold on
% pdfplot(fl,'g')
% pdfplot(f),pdfplot(fk,'r'),hold off
%
% See also  thvpdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway. 
  
% History
% By pab 10.01.2004
%error(nargchk(3,6,nargin))
narginchk(3,6)
if (nargin < 6||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end
if (nargin < 5||isempty(gam))
   gam = getjonswappeakedness(Hm0,Tp);
end
displayWarning = 1;
if displayWarning
  if any(any(Tp>5*sqrt(Hm0) | Tp<3.6*sqrt(Hm0)))
    disp('Warning: Hm0,Tp is outside the JONSWAP range')
    disp('The validity of the parameters returned are questionable')
  end
end

[V,H] = meshgrid(Vcf,Hd);

f = createpdf(2);
[f.f,Hrms,Vrms,varargout{1:nargout-1}]  = jhvnlpdf(H,V,Hm0,Tp,gam,normalizedInput);
f.x = {Vcf(:),Hd(:)};
 
if (normalizedInput)
  f.labx = {'Vcf', 'Hd'};
  f.norm = 1;
else
  f.norm = 0;
  f.labx = {'Vcf [m/s]', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Vcf) in time';
f.note = ['Jonswap Hm0=' num2str(Hm0) ' Tp = ' num2str(Tp) ' Gamma = ' num2str(gam)];
[f.cl,f.pl] = qlevels(f.f);
return