function [A0,B0,C0]= thwparfun(Hm0,Tp,dim)
%THWPARFUN Wave height, Hd, distribution parameters for Torsethaugen spectra.
%
% CALL [a b c] = thwparfun(Hm0,Tp,dim)
%
% Hm0 = significant wave height [m].
% Tp  = peak period [s]
% dim = 'time'  : Hd distribution parameters in time (default)
%
%  THWPARFUN returns the truncated Weibull distribution parameters which
%  approximates the marginal PDF of Hd/Hrms, i.e.,
%  zero-downcrossing wave height, for a Gaussian process with a 
%  Torsethaugen spectral density (torsethaugen).  
%
%  The empirical parameters of the model is fitted by
%  least squares to simulated Hd data for 600 classes of Hm0 and Tp.
%  Between 50000 and 150000 zero-downcrossing waves were simulated for
%  each class of Hm0 and Tp for DIM=='time'.
%  THWPARFUN is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
% 
%  Example:
%  Hm0 = 6;Tp = 8;Hrms = Hm0/sqrt(2);
%  [a b c] = thwparfun(Hm0,Tp);
%  h = linspace(0,4*Hrms)'; 
%  F = cdfweibmod(h/Hrms,a,b,c);
%  f = pdfweibmod(h/Hrms,a,b,c)/Hrms;
%  dt = 0.4; w = linspace(0,2*pi/dt,256)';
%  S = torsethaugen(w,[Hm0 Tp]);
%  xs = spec2sdat(S,20000,dt); rate=8; method=1;
%  [S,H] = dat2steep(xs,rate,method);
%  fk = kdebin(H,{'kernel','epan','L2',.5,'inc',128});
%  subplot(2,1,1)
%  plotedf(H,[h,F],1)
%  subplot(2,1,2)
%  plot(h,f), hold on, pdfplot(fk,'r'), hold off
% 
%  See also  thpdf 

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.
  
% History:
% revised pab 7-Jan-2004  
% by pab 29.11.2002

%error(nargchk(2,3,nargin))
narginchk(2,3)
if nargin<3||isempty(dim), dim = 'time';end

displayWarning = 0;
if displayWarning
  if  any(Hm0<0.5 | 12<Hm0)
    disp('Warning: Hm0 is outside the valid range')
    disp('The validity of the parameters returned are questionable')
  end
  if any(Tp<3.5 | 20<Tp)
    disp('Warning: Tp is outside the valid range')
    disp('The validity of the parameters returned are questionable')
  end

  if any(Hm0 > (Tp-2)*12/11)
    disp('Warning: Hm0 is too large compared to Tp!')
    disp('The validity of the parameters returned are questionable')
  end
end
global THSSPARW
global THWPAR
if strcmpi(dim(1),'t'), % Waveheight distribution in time
 
  if isempty(THWPAR)
    %THWPAR = load('thwpar.mat');
    THWPAR = load('thwnlpar20-Jul-2004.mat');
  end
  % Truncated Weibull  distribution parameters as a function of Tp, Hm0 
  A00 = THWPAR.A00s;
  B00 = THWPAR.B00s;
  C00 = THWPAR.C00s;

  Tpp  = THWPAR.Tp;
  Hm00 = THWPAR.Hm0;
  [E1, H1] = meshgrid(Tpp,Hm00);
  method = '*cubic';
  A0 = interp2(E1,H1,A00,Tp,Hm0,method);
  B0 = interp2(E1,H1,B00,Tp,Hm0,method);
  C0 = interp2(E1,H1,C00,Tp,Hm0,method);
else 
  % % Waveheight distribution in space not available
  %  disp('Waveheight distribution parameters in space not available')
  %  disp('for truncated weibull distribution.')
  %  disp('Corresponding distribution parameters is available for the gamma')
  %  disp('distribution in thgparfun!')
  
 
  if isempty(THSSPARW)
    THSSPARW = load('thsspar27-Jul-2004.mat');
  end
   % truncated Weibul distribution parameters as a function of Tp, Hm0 
   A00 = THSSPARW.A00s;
   B00 = THSSPARW.B00s;
   C00 = THSSPARW.C00s;
   Tpp  = THSSPARW.Tp;
   Hm00 = THSSPARW.Hm0;
   [E1, H1] = meshgrid(Tpp,Hm00);
   method = '*cubic';
   
   %A0 = repmat('nan',Hm0);  B0 = A0;  C0 = A0;
  
   A0 = interp2(E1,H1,A00,Tp,Hm0,method);
   B0 = interp2(E1,H1,B00,Tp,Hm0,method);
   C0 = interp2(E1,H1,C00,Tp,Hm0,method);
end
