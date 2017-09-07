function [f,varargout] = thsnlpdf2(Hd,Scf,Hm0,Tp,normalizedInput)
%THSNLPDF2 Joint (Scf,Hd) PDF for nonlinear waves with Torsethaugen spectra.
%
%  CALL: f = thsnlpdf2(Hd,Scf,Hm0,Tp)
% 
%   f   = pdf struct evaluated at meshgrid(Scf,Hd)
%   Hd  = zero down crossing wave height (vector)
%   Scf = crest front steepness (vector) 
%   Hm0 = significant wave height [m]
%   Tp  = Spectral peak period    [s]
%
% THSNLPDF2 approximates the joint distribution of (Scf, Hd), i.e., crest
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for 2nd order
% nonlinear waves with a
% Torsethaugen spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Scf,Hd) data for 600 classes of
% Hm0 and Tp. Between 40000 and 200000 zero-downcrossing waves were
% simulated for each class of Hm0 and Tp.
% THSNLPDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2)); 
% s = linspace(0,6*1.25*Hm0/Tp^2);
% f = thsnlpdf2(h,s,Hm0,Tp);
% w = linspace(0,40,5*1024+1).';
% S = torsethaugen(w,[Hm0 Tp]);
% dt = 0.3;  
% x = spec2nlsdat(S,80000,.2); rate = 8;
% [si,hi] = dat2steep(x,rate,2);
% fk = kdebin([si,hi],{'L2',.5,'inc',128});
%  fk.title = f.title; fk.labx = f.labx; 
% plot(si,hi,'.'), hold on
% pdfplot(f),pdfplot(fk,'r'),hold off
%
% See also  thsspdf

  
% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
% History
% revised pab 09.08.2003
% changed input  
% validated 20.11.2002
% By pab 20.12.2000

%error(nargchk(4,5,nargin))
narginchk(4,5)
if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end


[V,H] = meshgrid(Scf,Hd);

f = createpdf(2);
[f.f,Hrms,Vrms,varargout{1:nargout-1}]  = thsnlpdf(H,V,Hm0,Tp,normalizedInput);

 f.x = {Scf(:),Hd(:)};
 
if (normalizedInput)
  f.labx={'Scf', 'Hd'};
  f.norm = 1;
else
  f.norm=0;
  f.labx={'Scf', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Scf) in time (non-linear)';
f.note = ['Torsethaugen Hm0=' num2str(Hm0) ' Tp = ' num2str(Tp)];
[f.cl,f.pl] = qlevels(f.f);

return 
