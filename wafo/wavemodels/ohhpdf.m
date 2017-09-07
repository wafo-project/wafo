function f = ohhpdf(h,Hm0,def,dim)
%OHHPDF Marginal wave height, Hd, PDF for Bimodal Ochi-Hubble spectra. 
%
%  CALL: f = ohhpdf(h,Hm0,def,dim)
% 
%  f   = pdf evaluated at h.
%  h   = vectors of evaluation points.
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
% dim = 'time'  : Hd distribution in time (default)
%       'space' : Hd distribution in space
%
% OHHPDF approximates the marginal PDF of Hd, i.e.,
% zero-downcrossing wave height, for a Gaussian process with a Bimodal
% Ochi-Hubble spectral density (ochihubble). The empirical parameters of
% the model is fitted by least squares to simulated Hd data for 24
% classes of Hm0. Between 50000 and 150000 zero-downcrossing waves were
% simulated for each class of Hm0 in DIM=='time'.
% Between 50000 and 300000 zero-downcrossing waves were
% simulated for each class of Hm0 for DIM=='space'.
% OHHPDF is restricted to the following range for Hm0: 
%  0 < Hm0 [m] < 12,  1 <= def < 11, 
%
% Example:
% Hm0 = 6;def = 8; dim = 'time';
% h = linspace(0,4*Hm0/sqrt(2))'; 
% f = ohhpdf(h,Hm0,def,dim);
% plot(h,f)
% dt = 0.4; w = linspace(0,2*pi/dt,256)';
% xs = spec2sdat(ochihubble(w,[Hm0, def]),26000); rate=8; method=1;
% [S,H] = dat2steep(xs,rate,method);
% fk = kdebin(H,{'L2',.5,'inc',128}); 
% hold on, pdfplot(fk,'g'), hold off
%
% See also  ohhcdf,ohhvpdf

  
% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
% History
% revised pab Jan 2004
% By pab 20.01.2001


%error(nargchk(2,4,nargin))
narginchk(2,4)
if nargin<4||isempty(dim), dim  = 'time';end 
if nargin<3||isempty(def), def  = 1;end 

if Hm0>12|| Hm0<=0 
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end

if def>11||def<1 
  Warning('DEF is outside the valid range')
  def = mod(def-1,11)+1;
end

Hrms = Hm0/sqrt(2);

[A0 B0 C0] = ohhgparfun(Hm0,def,dim);
f = pdfgengam(h/Hrms,A0,B0,C0)/Hrms;
return

