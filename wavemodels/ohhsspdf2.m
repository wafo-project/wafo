function [f,varargout] = ohhsspdf2(Hd,Scf,Hm0,def,normalizedInput)
%OHHSSPDF2 Joint (Scf,Hd) PDF for linear waves in space with Ochi-Hubble spectra. 
%
%  CALL: f = ohhsspdf2(Hd,Scf,Hm0,Tp)
% 
%  f   = pdf struct evaluated at meshgrid(Scf,Hd).
%  Hd  = zero down crossing wave height
%  Scf = crest front steepness
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
%
% OHHSSPDF2 approximates the joint distribution of (Scf, Hd), i.e., crest
% front steepness (Ac/Lcf) and wave height in space, for a Gaussian
% process with a bimodal Ochi-Hubble spectral density. The empirical
% parameters of the model is fitted by least squares to simulated
% (Scf,Hd) data for 24 classes of Hm0.
% Between 50000 and 300000 zero-downcrossing waves were simulated for
% each class of Hm0.
% OHHSSPDF2 is restricted to the following range for Hm0: 
%  0.5 < Hm0 [m] < 12
%
% Example:
% Hm0 = 6;Tp = 8;def= 2;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,Hm0/Tp)';
% f = ohhsspdf2(h,v,Hm0,def);
% w = linspace(0,10,2*1024+1).';
% S = ochihubble(w,[Hm0 def]);
% Sk = spec2spec(specinterp(S,.55),'k1d');
% dk = 1;
% x = spec2sdat(Sk,80000,dk); rate = 8;
% [vi,hi] = dat2steep(x,rate,1);
% fk = kdebin([vi,hi],{'L2',.5,'inc',128});
% fk.title = f.title; fk.labx = f.labx; 
% plot(vi,hi,'.'), hold on
% pdfplot(f),
% pdfplot(fk,'r'), hold off
%
% See also  ohhspdf, thvpdf

  
% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
% History
% revised pab 09.09.2003
% By pab 20.12.2000


error(nargchk(4,5,nargin))

if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(def)),  def  = 1;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end


[V,H] = meshgrid(Scf,Hd);

f = createpdf(2);
[f.f,Hrms,Vrms,varargout{1:nargout-1}]  = ohhsspdf(H,V,Hm0,def,normalizedInput);

 f.x = {Scf(:),Hd(:)};
 
if (normalizedInput)
  f.labx={'Scf', 'Hd'};
  f.norm = 1;
else
  f.norm=0;
  f.labx={'Scf', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Scf) in space';
f.note = ['ohhspec2 Hm0=' num2str(Hm0) ' def = ' num2str(def)];
[f.cl,f.pl] = qlevels(f.f);

return
