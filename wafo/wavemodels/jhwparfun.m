function [A0,B0,C0]= jhwparfun(Hm0,Tp,gam,dim)
%JHWPARFUN Wave height, Hd, distribution parameters for Jonswap spectrum. 
%
% CALL [a b c] = jhwparfun(Hm0,Tp,gamma,dim)
%
% Hm0   = significant wave height [m].
% Tp    = peak period [s]
% gamma = Peakedness parameter of the JONSWAP spectrum
% dim   = 'time'  : Hd distribution parameters in time (default)
%
%  JHWPARFUN returns the truncated Weibull distribution parameters which
%  approximates the marginal PDF of Hd/Hrms, i.e.,
%  zero-downcrossing wave height, for a Gaussian process with a 
%  Jonswap spectral density.  
%
%  The empirical parameters of the model is
%  fitted by least squares to simulated Hd data for 13 classes of
%  GAMMA between 1 and 7.
%  About 50000 zero-downcrossing waves were simulated for
%  each class of GAMMA for DIM=='time'.
%  JHWPARFUN is restricted to the following range for GAMMA and Tp: 
%  1 <= GAMMA <= 7, and 3.6*sqrt(Hm0) < Tp < 5*sqrt(Hm0)
% 
%  Example:
%  Hm0 = 6;Tp = 9;Hrms = Hm0/sqrt(2);
%  gam = getjonswappeakedness(Hm0,Tp);
%  [a b c] = jhwparfun(Hm0,Tp,gam);
%  h = linspace(0,4*Hrms)'; 
%  F = cdfweibmod(h/Hrms,a,b,c);
%  f = pdfweibmod(h/Hrms,a,b,c)/Hrms;
%  dt = 0.4; w = linspace(0,2*pi/dt,256)';
%  S = jonswap(w,[Hm0 Tp,gam]);
%  xs = spec2sdat(S,80000,dt); rate=8; method=1;
%  [S,H] = dat2steep(xs,rate,method);
%  fk = kdebin(H,{'L2',.5,'inc',128});
%  subplot(2,1,1)
%  plotedf(H,[h,F],1)
%  subplot(2,1,2)
%  plot(h,f), hold on, pdfplot(fk,'r'), hold off
% 
%  See also  jhvpdf 

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.
  
% History:
% revised pab 4 jan 2006
%  updated example in help header
% revised pab 10 jan 2004  
% by pab 29.11.2002

%error(nargchk(2,4,nargin))
narginchk(2,4)
if nargin<3||isempty(gam),
  gam = getjonswappeakedness(Hm0,Tp);
end
if nargin<4||isempty(dim), dim = 'time';end

displayWarning = 0;
if displayWarning
  if any(any(Tp>5*sqrt(Hm0) | Tp<3.6*sqrt(Hm0)))
    disp('Warning: Hm0,Tp is outside the JONSWAP range')
    disp('The validity of the parameters returned are questionable')
  end
  if any(any(gam>7|gam<1))
    disp('Warning: gamma is outside the valid range')
    disp('The validity of the parameters returned are questionable')
  end
end

if strncmpi(dim,'t',1)    
  % LS fit to data
  % best fit to jonswap for gamma = 1:.5:7
  A0 =-0.00768035933493*gam +  1.12723230522407;
  B0 =-0.02261969635571*gam +  2.31625911055348;
  C0 = -0.00808703761189*gam +  0.14438532196018; 
else % not implemented yet
  
  A0 = [];
  B0 = [];
  C0 = [];
end

return

% A0 = -0.01243795213128.*gam + 1.08025514722235;
% B0 = -0.03043834819688.*gam + 2.27161821064622;
% C0 = -0.01342959276544.*gam + 0.10353423379696;
