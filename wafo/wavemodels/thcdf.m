function [f] = thcdf(h,Hm0,Tp,dim)
%THCDF Marginal wave height, Hd, cdf for Torsethaugen spectra. 
%
%  CALL: F = thcdf(h,Hm0,Tp)
% 
%  F   = cdf evaluated at h.
%  h   = vector of evaluation points.
%  Hm0 = significant wave height [m].
%  Tp  = Spectral peak period    [s].
%  dim = 'time'  : Hd distribution in time (default)
%        'space' : Hd distribution in space
%
% THCDF approximates the marginal cumulative distribution of Hd, i.e.,
% zero-downcrossing wave height, for a Gaussian process with a Torsethaugen
% spectral density. The empirical parameters of the model is fitted by
% least squares to simulated Hd data for 600 classes of Hm0 and
% Tp. Between 50000 and 150000 zero-downcrossing waves were simulated for
% each class of Hm0 and Tp.
% THCDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% F = thcdf(h,Hm0,Tp);
% dt = 0.4; w = linspace(0,2*pi/dt,256)';
% S = torsethaugen(w,[Hm0 Tp]);
% xs = spec2sdat(S,20000,dt); rate=8; method=1;
% [S,H] = dat2steep(xs,rate,method);
% plotedf(H,[h, F],'g')
%
% See also  thpdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.
  
% History
% Revised pab jan2004  
% By pab 20.12.2000


error(nargchk(3,4,nargin))
if nargin<4||isempty(dim),
  dim = 'time';% dim='time'->wtweibpdf, dim='space'->wggampdf
end

if any( Hm0>12| Hm0>(Tp-2)*12/11)
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end
if any(Tp>20|Tp<3) 
  disp('Warning: Tp is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end
Hrms = Hm0/sqrt(2);
[a b c] = thwparfun(Hm0,Tp,dim);
f = cdfweibmod(h/Hrms,a,b,c);
return

%old call kept just in case
% if strncmpi(dim,'t',1)    
%   [a b c] = thwparfun(Hm0,Tp);
%   f = cdfweibmod(h/Hrms,a,b,c);
% else 
%   [a b c] = thgparfun(Hm0,Tp,dim);
%   f = cdfgengam(h/Hrms,a,b,c);
% end
% % return

