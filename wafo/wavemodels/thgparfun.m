function [A0,B0,C0]= thgparfun(Hm0,Tp,dim)
%THGPARFUN Wave height, Hd, distribution parameters for Torsethaugen spectra. 
%
% CALL [a b c] = thgparfun(Hm0,Tp,dim)
%
% Hm0 = significant wave height [m].
% Tp  = peak period [s]
% dim = 'time'  : Hd distribution parameters in time (default)
%       'space' : Hd distribution parameters in space
%
%  THGPARFUN returns the Generalized gamma distribution parameters which
%  approximates the marginal PDF of Hd/Hrms, i.e.,
%  zero-downcrossing wave height, for a Gaussian process with a 
%  Torsethaugen spectral density (torsethaugen).  
%
%  The empirical parameters of the model is fitted by
%  least squares to simulated Hd data for 600 classes of Hm0 and Tp.
%  Between 50000 and 150000 zero-downcrossing waves were simulated for
%  each class of Hm0 and Tp for DIM=='time'.
%  Between 100000 and 1000000 zero-downcrossing waves were
%  simulated for each class of Hm0 and Tp for DIM=='space'.
%  THGPARFUN is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
% 
%  Example:
%  Hm0 = 6;Tp = 8;Hrms = Hm0/sqrt(2);
%  [a b c] = thgparfun(Hm0,Tp);
%  h = linspace(0,4*Hrms)'; 
%  F = cdfgengam(h/Hrms,a,b,c);
%  f = pdfgengam(h/Hrms,a,b,c)/Hrms;
%  dt = 0.4; w = linspace(0,2*pi/dt,256)';
%  S = torsethaugen(w,[Hm0 Tp]);
%  xs = spec2sdat(S,20000,dt); rate=8; method=1;
%  [S,H] = dat2steep(xs,rate,method);
%  subplot(2,1,1)
%  plotedf(H,[h,F],1)
%  fk = kdebin(H,{'kernel','epan','L2',.5,'inc',128});
%  subplot(2,1,2)
%  plot(h,f), hold on, pdfplot(fk,'r'), hold off
% 
%  See also  thpdf 


% History:
% by pab 29.11.2002

%error(nargchk(2,3,nargin))
narginchk(2,3)
persistent THGPAR  THSSPAR

if nargin<3||isempty(dim), dim = 'time';end

if  any(Hm0<=0.5 | 12<Hm0)
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

pardef = 1;
switch pardef
  case 1,
    if strncmpi(dim,'t',1), % Waveheight distribution in time
      disp('Note: truncated weibull distribution gives better fit')
      disp('than this, see thwparfun for details!')
      
      
      if isempty(THGPAR)
        THGPAR = load('thgpar.mat');
      end
      % Generalized Gamma  distribution parameters as a function of Tp, Hm0 
      A00 = THGPAR.A00s;
      B00 = THGPAR.B00s;
      C00 = THGPAR.C00s;

      Tpp  = THGPAR.Tp;
      Hm00 = THGPAR.Hm0;
    else % Waveheight distribution in space
      %global THSSPAR
      if isempty(THSSPAR)
        THSSPAR = load('thsspar.mat');
      end
      % Generalized Gamma  distribution parameters as a function of Tp, Hm0 
      A00 = THSSPAR.A00s;
      B00 = THSSPAR.B00s;
      C00 = THSSPAR.C00s;

      Tpp  = THSSPAR.Tp;
      Hm00 = THSSPAR.Hm0;     
    end
    [E1, H1] = meshgrid(Tpp,Hm00);
    method = '*cubic';
    A0 = interp2(E1,H1,A00,Tp,Hm0,method);
    B0 = interp2(E1,H1,B00,Tp,Hm0,method);
    C0 = interp2(E1,H1,C00,Tp,Hm0,method);
end
