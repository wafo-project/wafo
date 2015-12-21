function [f,varargout] = thvpdf2(Hd,Vcf,Hm0,Tp,normalizedInput)
%THVPDF2 Joint (Vcf,Hd) PDF for linear waves with Torsethaugen spectra. 
%
%  CALL: f = thvpdf2(Hd,Vcf,Hm0,Tp)
% 
%  f   = pdf structure evaluated at meshgrid(Vcf,Hd).
%  Hd  = zero down crossing wave height
%  Vcf = crest front velocity
%  Hm0 = significant wave height [m].
%  Tp  = Spectral peak period    [s].
%
% THVPDF2 approximates the joint PDF of (Vcf, Hd), i.e., crest
% front velocity (Ac/Tcf) and wave height, for a Gaussian process with a
% Torsethaugen spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 600 classes of
% Hm0 and Tp. Between 50000 and 150000 zero-downcrossing waves were
% simulated for each class of Hm0 and Tp.
% THVPDF2 is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;
%  Tp = 8;
%  h  = linspace(0,4*Hm0/sqrt(2))'; 
%  v  = linspace(0,4*2*Hm0/Tp)';
%  f  = thvpdf2(h,v,Hm0,Tp);
%  w  = linspace(0,40,5*1024+1).';
%  S  = torsethaugen(w,[Hm0 Tp]);
%  dt = 0.3;
%  x  = spec2sdat(S,80000,.2); rate = 8;
%  [vi,hi] = dat2steep(x,rate,1);
%  fk = kdebin([vi,hi],{'L2',.5,'inc',128});
%  fk.title = f.title; fk.labx = f.labx; 
%  plot(vi,hi,'.'), hold on
%  pdfplot(f),pdfplot(fk,'r'),hold off
%
% See also  thspdf, thsspdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
% History
% revised pab 10.08.2003  
% -revised pab 28.11.2002
%   extended the example
% By pab 20.12.2000

error(nargchk(3,5,nargin))

if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end

[V,H] = meshgrid(Vcf,Hd);

f = createpdf(2);
[f.f,Hrms,Vrms,varargout{1:nargout-1}]  = thvpdf(H,V,Hm0,Tp,normalizedInput);
f.x = {Vcf(:),Hd(:)};
 
if (normalizedInput)
  f.labx = {'Vcf', 'Hd'};
  f.norm = 1;
else
  f.norm = 0;
  f.labx = {'Vcf [m/s]', 'Hd [m]'};
end
f.title = 'Joint distribution of (Hd,Vcf) in time';
f.note = ['Torsethaugen Hm0=' num2str(Hm0) ' Tp = ' num2str(Tp)];
[f.cl,f.pl] = qlevels(f.f);
return
