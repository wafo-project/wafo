function [f,Hrms,Srms] = jhspdf2(Hd,Scf,Hm0,Tp,gam,normalizedInput)
%JHSPDF2 Joint (Scf,Hd) PDF for linear waves with a JONSWAP spectrum.
%
%  CALL: f = jhspdf2(Hd,Scf,Hm0,Tp,gamma)
% 
%  f     = pdf struct evaluated at meshgrid(Scf,Hd)
%  Hd    = zero down crossing wave height
%  Scf   = crest front steepness
%  Hm0   = significant wave height
%  Tp    = Spectral peak period 
%  Gamma = Peakedness parameter of the JONSWAP spectrum
%
% JHSPDF2 approximates the joint distribution of (Scf, Hd), i.e., crest
% front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for a Gaussian
% process with a JONSWAP spectral density. The empirical parameters of
% the model is fitted by least squares to simulated (Scf,Hd) data for 13
% classes of GAMMA between 1 and 7. Between 47000 and 55000
% zero-downcrossing waves were simulated for each class.
% JHSPDF2 is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%
% Example:
% Hm0 = 6;Tp = 9; gam=3.5
% h = linspace(0,4*Hm0/sqrt(2))'; 
% s = linspace(0,6*1.25*Hm0/Tp^2)';
% f = jhspdf2(h,s,Hm0,Tp,gam);
% w = linspace(0,40,5*1024+1).';
% S = jonswap(w,[Hm0, Tp, gam]);
% dt = .3;
% x = spec2sdat(S,80000,dt); rate = 4;
% [si,hi] = dat2steep(x,rate,2);
% fk = kdebin([si,hi],{'kernel','epan','L2',.5,'inc',128}); 
% fk.title = f.title; fk.labx = f.labx;  
% plot(si,hi,'.'), hold on
% pdfplot(f),pdfplot(fk,'r'),hold off
%
% See also  jhspdf

  
% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway. 
  
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

[V,H] = meshgrid(Scf,Hd);

f = createpdf(2);
[f.f,Hrms,Srms,varargout{1:nargout-1}]  = jhspdf(H,V,Hm0,Tp,gam,normalizedInput);
f.x = {Scf(:),Hd(:)};
 
if (normalizedInput)
  f.labx = {'Scf', 'Hd'};
  f.norm = 1;
else
  f.norm = 0;
  f.labx = {'Scf [m/s]', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Scf) in time';
f.note = ['Jonswap Hm0=' num2str(Hm0) ' Tp = ' num2str(Tp) ' Gamma = ' num2str(gam)];
[f.cl,f.pl] = qlevels(f.f);
return
