function f = ohhcdf(h,Hm0,def,dim)
%OHHCDF Marginal wave height, Hd, CDF for Ochi-Hubble spectra.
%
%  CALL: f = ohhcdf(h,Hm0,def,dim)
% 
%  f   = cdf evaluated at h.
%  h   = vectors of evaluation points.
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
% dim = 'time'  : Hd distribution in time (default)
%       'space' : Hd distribution in space
%
% OHHCDF approximates the marginal PDF of Hd, i.e.,
% zero-downcrossing wave height, for a Gaussian process with a Bimodal
% Ochi-Hubble spectral density (ochihubble). The empirical parameters of
% the model is fitted by least squares to simulated Hd data for 24
% classes of Hm0. Between 50000 and 150000 zero-downcrossing waves were
% simulated for each class of Hm0 in DIM=='time'.
% Between 50000 and 300000 zero-downcrossing waves were
% simulated for each class of Hm0 for DIM=='space'.
% OHHCDF is restricted to the following range for Hm0: 
%  0 < Hm0 [m] < 12,  1 <= def < 11, 
%
% Example:
% Hm0 = 6;def = 8;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% f = ohhcdf(h,Hm0,def);
% plot(h,f)
% dt = 0.4; w = linspace(0,2*pi/dt,256)';
% xs = spec2sdat(ochihubble(w,[Hm0, def]),6000); rate=8; method=1;
% [S,H] = dat2steep(xs,rate,method);
% plotedf(H,[h f],'g')
%
% See also  ohhpdf, ohhvpdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
% Adapted to  cssmooth  by GL Feb 2011
% History
% revised pab jan2004  
% By pab 20.01.2001


%error(nargchk(2,4,nargin))
narginchk(2,4)
if nargin<4||isempty(dim), dim  = 'time';end 
if nargin<3||isempty(def), def  = 1;end 

if any(Hm0>12| Hm0<=0) 
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end

if def>11||def<1 
  warning('DEF is outside the valid range')
  def = mod(def-1,11)+1;
end

Hrms = Hm0/sqrt(2);
[A0, B0, C0] = ohhgparfun(Hm0,def,dim);
f    = cdfgengam(h/Hrms,A0,B0,C0);
return
%old calls
% pardef = 7;
% switch pardef
%   case 1
%      w    = linspace(0,100,16*1024+1).'; % original spacing
%      S = ochihubble(w,[Hm0,def]);
%      R  = spec2cov(S);
%      %    A0 = sqrt((1-min(R.R)/R.R(1))/2);% Naess (1985)
%      A0 = sqrt((1-min(R.R)/R.R(1))/2)+0.03;% Modified approach broadbanded time
%      %    A0 = sqrt((1-min(R.R)/R.R(1))/2)+0.1;% Modified approach broadbanded space                                 
					  
%     B0 = 2;
%     C0 = 0;
%   case 7,
%     global OHHWPAR
%     if isempty(OHHWPAR)
%       OHHWPAR = load('thwpar.mat');
%     end
%     % Truncated Weibull  distribution parameters as a function of Tp, Hm0 
%     A00 = OHHWPAR.A00s;
%     B00 = OHHWPAR.B00s;
%     C00 = OHHWPAR.C00s;

%     Hm00 = OHHWPAR.Hm0;
%     if 1,
%       method = '*cubic';
%       A0 = interp1(Hm00,A00(:,def),Hm0,method);
%       B0 = interp1(Hm00,B00(:,def),Hm0,method);
%       C0 = interp1(Hm00,C00(:,def),Hm0,method);
%     else
%       A0 = cssmooth(Hm00,A00(:,def),1,Hm0);
%       B0 = cssmooth(Hm00,B00(:,def),1,Hm0);
%       C0 = cssmooth(Hm00,C00(:,def),1,Hm0);
%     end
 
% end

%Hrms = Hm0/sqrt(2);
%f.f    = cdfweibmod(h/Hrms,A0,B0,C0)/Hrms;


