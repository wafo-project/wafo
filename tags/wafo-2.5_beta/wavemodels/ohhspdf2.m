function [f,varargout] = ohhspdf2(Hd,Scf,Hm0,def,normalizedInput)
%OHHSPDF2 Joint (Scf,Hd) PDF for linear waves with Ochi-Hubble spectra. 
%
%  CALL: f = ohhspdf2(Hd,Scf,Hm0,def)
% 
%  f   = pdf struct evaluated at meshgrid(Scf,Hd)
%  Hd  = zero down crossing wave height
%  Scf = crest front steepness
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
%
% OHHSPDF2 approximates the joint distribution of (Scf, Hd) in time, 
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height,
%  for a Gaussian process with a bimodal Ochi-Hubble spectral density
% (ochihubble). The empirical
% parameters of the model is fitted by least squares to simulated
% (Scf,Hd) data for 24 classes of Hm0. Between 50000 and 150000
% zero-downcrossing waves were simulated for each class of Hm0.
% OHHSPDF is restricted to the following range for Hm0: 
% 0.5 < Hm0 [m] < 12
% The size of f is the common size of the input arguments, Hd, Scf and
% Hm0.  
%
% Example:
% Hm0 = 6;Tp = 8;def= 2;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,6*1.25*Hm0/Tp^2)';
% f = ohhspdf2(h,v,Hm0,def);
% w = linspace(0,40,5*1024+1).';
% S = ochihubble(w,[Hm0 def]);
% dt = 0.3;
% x = spec2sdat(S,80000,.2); rate = 8;
% [si,hi] = dat2steep(x,rate,2);
% fk = kdebin([si,hi],{'L2',.5,'inc',128});
%  fk.title = f.title; fk.labx = f.labx; 
% plot(si,hi,'.'), hold on
% pdfplot(f),pdfplot(fk,'r'),hold off  
%
% See also  ohhpdf, thspdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
  
% History
% revised pab 09.09.2003
% By pab 06.02.2001



error(nargchk(4,5,nargin))

if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(def)),  def  = 1;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end


[V,H] = meshgrid(Scf,Hd);

f = createpdf(2);
[f.f,Hrms,Vrms,varargout{1:nargout-1}]  = ohhspdf(H,V,Hm0,def,normalizedInput);

 f.x = {Scf(:),Hd(:)};
 
if (normalizedInput)
  f.labx={'Scf', 'Hd'};
  f.norm = 1;
else
  f.norm=0;
  f.labx={'Scf', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Scf) in time';
f.note = ['ohhspec2 Hm0=' num2str(Hm0) ' def = ' num2str(def)];
[f.cl,f.pl] = qlevels(f.f);

return
