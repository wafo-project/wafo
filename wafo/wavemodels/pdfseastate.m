function p = pdfseastate(t,h,Pt1,Ph1, def,condon)
%PDFSEASTATE Long Term Joint Sea State PDF of Hm0 and wave period 
%
% CALL: f = pdfseastate(T,Hs,Pt,Ph,def)
%
%    f      = density 
%    T      = peak period, Tp  or mean zerodowncrossing period, Tz (sec)
%    Hs     = significant wave height (meter)
%    Pt     = [A1 A2 A3 B1 B2 B3 B4] parameters for the conditional
%             distribution of T given Hs 
%             (default [1.59 0.42 2 0.005 0.09 0.13 1.34])
%    Ph     = [C1 C2 C3 Hc M1 S1], parameters for the marginal
%             distribution for Hs (default [2.82 1.547 0 3.27  0.836 0.376])
%    def    = defines the parametrization used (default 0)
%
%   The probability distribution is given as the product:
%
%    f(T,Hs)=f(T|Hs)*f(Hs)
%   where
%
%       f(T|Hs) = pdflognorm(T,my(Hs),S(Hz)) with mean and variance
%        my(Hs)  = A1 + A2*log(Hs+A3)          if mod(def,2)==0
%                  A1 + A2*Hs^A3               if mod(def,2)==1
%        S(Hs)   = B1+B2*exp(-B3*Hs.^B4)
%
%    For def <=1: 
%       f(Hs)    = pdflognorm(Hs,M1,S1)    for Hs <= Hc 
%                = pdfweib(Hs-C3,C1,C2)    for Hs >  Hc  
%    For def >1:   
%       f(Hs)    = pdflognorm(Hs,M1,S1)    for Hs <= Hc 
%                = pdfgengam(Hs,C1,C2,C3)  for Hs >  Hc  (Suggested by Ochi)
%
%    The default values for T and Hs are suitable for the Northern North
%    Sea for peak period, Tp and significant waveheight, Hs.
%    With a suitable change of parameters for Pt this model fit mean wave 
%    period, Tz, also reasonably well.
%    
%   The size of f is the common size of T and Hs
%
%   NOTE: - by specifying nan's  in the vectors Pt or Ph default values
%           will be used.
%         - if length(Pt) or length(Ph) is shorter than the parameters
%           needed then the default values are used for the parameters
%           not specified. 
%         - For tables of fitted parameter values to specific sites see
%           the end of this file (type pdfseastate.m)
%
% Example: % Set  C1 = 2.73,  Hc = 3.95  and the rest to their default values
%   Ph = [2.73, NaN,NaN, 3.95]; x = linspace(0, 15);
%   [T,H] = meshgrid(x);
%   f = pdfseastate(T,H,[],Ph);
%   contour(x,x,f)
%   
% See also  pdfweib, pdflognorm, pdfgengam

%   References:
%   Haver, S (1980) 
%  'Analysis of uncertainties related to the stochastic modelling of
%    Ocean waves'
%   Ph.D. thesis, Norwegian Institute of Technology, NTH, Trondheim, Norway
%
%   Haver,S and Nyhus, K. A.  (1986) 
%   'A wave climate description for long term response calculation.'
%   In Proc. of OMAE'86, tokyo, Japan, pp. 27-34
%
%  Sagli, Gro (2000)
%  "Model uncertainty and simplified estimates of long term extremes of
%  hull girder loads in ships"
%  Ph.D. thesis, Norwegian University of Science and Technology, NTNU,
%  Trondheim, Norway, pp 45--47
%
%   Michel K. Ochi (1998),
%   "OCEAN WAVES, The stochastic approach",
%   OCEAN TECHNOLOGY series 6, Cambridge, pp 127.    
% (Generalized Gamma distribution for wave height)
%
%  Bitner-Gregersen, E.M. and Hagen, ุ. (2000)
% "Aspects of joint distribution for metocean Phenoma at the Norwegian
% Continental shelf."
% In Proc. of OMAE'2000, 

% tested on: matlab 5.2
% history:
% revised pab 2008
% renamed from ltwcpdf to pdfseastate
% revised pab 04.11.2000
%  -updated calls to new wstats functions
% by  Per A. Brodtkorb 13.05.2000





error(nargchk(2,6,nargin))

if nargin < 5||isempty(def),  def=0; end
if nargin < 6||isempty(condon),  condon=0; end

Ph=[2.82 1.547 0 3.27  0.836 0.376]; % default values
if nargin < 4||isempty(Ph),
 
else
  nh=length(Ph1); 
  ind=find(~isnan(Ph1(1:min(nh,6))));
  if any(ind) % replace default values with those from input data
    Ph(ind)=Ph1(ind);
  end
end

Pt=[1.59 0.42 2 0.005 0.09 0.13 1.34]; % default values
if nargin < 3||isempty(Pt),
else
  nt=length(Pt1); 
  ind=find(~isnan(Pt1(1:min(nt,7))));
  if any(ind) % replace default values with those from input data
    Pt(ind) = Pt1(ind);
  end
end

[icode,t,h] = iscomnsize(t,h);
if ~icode 
  error('Requires non-scalar arguments to match in size.');
end


% Model parameters for the marginal distribution for Hz:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% scale, shape and location parameter for weibull distribution
% or shape,shape and scale parameter for the generalized Gamma distribution
C1=Ph(1);C2=Ph(2);C3=Ph(3);
% Split value, mean and standard deviation for Lognormal distribution
Hc=Ph(4);myHz=Ph(5);sHz=Ph(6);

% Model parameters for the conditional distribution of Tp or Tz given Hz
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A1=Pt(1);A2=Pt(2);A3=Pt(3);
B1=Pt(4);B2=Pt(5);B3=Pt(6);B4=Pt(7);
% Mean and Variance of the lognormal distribution
if mod(abs(def),2)==0
  myTp = A1+A2*log(h+A3);
else
  myTp = A1+A2*h.^A3;
end
sTp =(B1+B2*exp(-B3*h.^B4));

p=zeros(size(h));
k0=find(t>0 & h>0);
if any(k0)
  if condon==5,
     p(k0) = cdflognorm(t(k0),myTp(k0),sTp(k0));
   elseif condon ==10
     p(k0) = 1;
   else
    p(k0) = pdflognorm(t(k0),myTp(k0),sTp(k0));
  end
end



if condon==0 || condon==5 || condon ==10
  k=find(0<h & h<=Hc );
  if any(k)
    p(k)=p(k).*pdflognorm(h(k),myHz,sHz);
  end
  k1=find(h>Hc);
  if any(k1)
    if def<=1,
      % NB! weibpdf must be modified to correspond to
      % pdf=x^(b-1)/a^b*exp(-(x/a)^b)
      p(k1)=p(k1).*pdfweib(h(k1)-C3,C1,C2);
    else
      p(k1)=p(k1).*pdfgengam(h(k1),C1,C2,C3); 
    end
  end
end

switch condon,
  case 0, % regular pdf is returned 
  case 1, %pdf conditioned on Hz ie. p(Tp|Hz) 
  case 3, % secret option  used by ltwcstat: returns Tp*p(Tp|Hz) 
    p = t.*p;
  case 4, % secret option  used by ltwcstat: returns Tp.^2*p(Tp|Hz) 
    p = t.^2.*p;
  case 5, % p(Hz)*P(Tp|Hz) is returned special case used by ltwccdf
  case 10 % p(Hz) is returned
   
  otherwise, error('unknown option')
end




% The following is taken from Sagli (2000):
% Parameters for the long term description of the wave climate (Tp|Hs):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area            Pt = [  A1     A2    A3     B1    B2    B3   B4]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   Northern North Sea:  [1.59,  0.42, 2.00, 0.005, 0.09, 0.13, 1.34]  0
%   Aasgard           :  [1.72,  0.34, 0.46, 0.005, 0.10, 0.29, 1.00]  1
%   Sleipner          :  [0.23,  1.69, 0.15, 0.005, 0.12, 0.40, 1.00]  1
%   Troms I           :  [1.35,  0.61, 0.34, 0.005, 0.12, 0.36, 1.00]  1
%   Statfjord/Brent   :  [0.40,  1.59, 0.15, 0.005, 0.12, 0.34, 1.00]  1
%   Ekofisk           :  [1.00,  0.90, 0.25, 0.005, 0.10, 0.44, 1.00]  1
%   Ekofisk Zone II   :  [0.03,  1.81, 0.15, 0.005, 0.16, 0.58, 1.00]  1
%   Ekofisk Zone III  :  [0.03,  1.81, 0.15, 0.005, 0.16, 0.58, 1.00]  1
%   Ekofisk Zone IV   :  [1.41,  0.38, 0.54, 0.010, 0.08, 0.41, 1.00]  1
%   Ekofisk Zone VI   :  [0.00,  1.72, 0.15, 0.010, 0.14, 1.00, 1.00]  1
%
%
% Parameters for the long term description of the wave climate (Hs):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area            Ph = [  C1     C2    C3    Hc    M1    S1  ]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   Northern North Sea:  [2.82,  1.55, 0.00, 3.27, 0.84, 0.61^2  ]  0
%   Aasgard           :  [2.73,  1.45, 0.00, 3.95, 0.84, 0.59^2  ]  1
%   Sleipner          :  [2.38,  1.44, 0.00, 3.00, 0.69, 0.60^2  ]  1
%   Troms I           :  [2.26,  1.36, 0.00, 4.25, 0.72, 0.56^2  ]  1
%   Statfjord/Brent   :  [2.82,  1.55, 0.00, 3.27, 0.84, 0.61^2  ]  1
%   Ekofisk           :  [2.08,  1.37, 0.00, 2.65, 0.52, 0.67^2  ]  1
%   Ekofisk Zone II   :  [2.12,  1.37, 0.00, 2.70, 0.54, 0.67^2  ]  1
%   Ekofisk Zone III  :  [2.05,  1.55, 0.00, 2.15, 0.50, 0.65^2  ]  1
%   Ekofisk Zone IV   :  [1.56,  1.46, 0.00, 2.00, 0.25, 0.62^2  ]  1
%   Ekofisk Zone VI   :  [1.18,  1.42, 0.00, 1.25, 0.07, 0.70^2  ]  1
%


% The following is taken from Bitner-Gregersen and Hagen (2000):
% Description of sites:
% V๘ring plateau:  (position 67, 27'N 5, 58'E water depth 1460m) NORWAVE
%                   buoy every 3rd hour for a period of 2.3 years 
% 

% Parameters for the long term description of the wave climate (Tp|Hs):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area           Pt = [  A1     A2      A3      B1     B2       B3     B4]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   V๘ring M         :  [1.718,  0.194, 0.642, -288.28, 288.52,0.279e-5 , 1]  1
%   V๘ring M/H       :  [1.670,  0.220, 0.617, - 93.13,  93.37,0.330e-4 , 1]  1
%   V๘ring H         :  [0.474,  1.438, 0.186,   0.047,   0.46,   0.476 , 1]  1
%   Haltenbanken     :  [1.921,  0.171, 0.665,   0.064,  0.327,   0.227 , 1]  1

% Parameters for the long term description of the wave climate (Tz|Hs)
% (mean zero-down crossing period given Hz)                           :
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area            Pt = [  A1     A2      A3     B1    B2     B3   B4]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   Statfjord&Gullfaks:  [1.790,  0.110, 0.759, 0.106, 0.361,0.969 , 1]  1
%   Statfjord Platform:  [1.771,  0.080, 0.902, 0.001, 0.239,0.169 , 1]  1
%   Ekofisk           :  [0.611,  0.902, 0.289, 0.071, 0.211,0.606 , 1]  1


% Parameters for the long term description of the wave climate (Hs):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area            Ph = [  C1     C2    C3    Hc  ]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   V๘ring M          :  [2.371,  1.384, 0.755, 0  ]  1
%   V๘ring M/H        :  [2.443,  1.405, 0.820, 0  ]  1
%   V๘ring H          :  [2.304,  1.383, 0.899, 0  ]  1
%   Haltenbanken      :  [2.154,  1.273, 0.763, 0  ]  1
%   Statfjord&Gullfaks:  [2.264,  1.398, 0.969, 0  ]  1
%   Statfjord Platform:  [2.502,  1.527, 0.657, 0  ]  1
%   Ekofisk           :  [1.597,  1.184, 0.852, 0  ]  1


% Fitted parameters by splitting the data into wind Sea and swell:

% Parameters for the long term description of the wave climate (Tp|Hs):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area            Pt = [  A1     A2    A3     B1    B2    B3   B4]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   V๘ring M         :  [-0.707, 2.421, 0.126, 0.076, 0.063, 0.162, 1]  1    (Wind Sea)
%   V๘ring M         :  [ 1.501, 0.882, 0,135, 0.126, 0.000, 0.000, 1]  1    (Swell)
%   Haltenbanken     :  [-4383.9, 4385.6, 0.948e-4, 0.045, 0.369,0.780,1]  1 (Wind Sea)
%   Haltenbanken     :  [ 2.257, 0.088, 0.746, 0.010, 0.229, 0.167, 1]  1    (Swell)

% Parameters for the long term description of the wave climate (Hs):
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Area            Ph = [  C1     C2    C3    Hc  ]  def    
%   จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
%   V๘ring M          :  [2.170,  1.382, 0.203, 0  ]  1  (Wind sea)
%   V๘ring M          :  [1.606,  1.252, 0.229, 0  ]  1  (Swell)
%   Haltenbanken      :  [1.473,  1.032, 0.241, 0  ]  1  (Wind sea)
%   Haltenbanken      :  [1.661,  1.117, 0.240, 0  ]  1  (Swell)




