function S1=pmspec(w1,sdata,plotflag)
%PMSPEC Calculates (and plots) a Pierson-Moskowitz spectral density.
% 
% CALL:  S = pmspec(w,data,plotflag); 
%        S = pmspec(wc,data,plotflag);
% 
%        S    = a struct containing the spectral density, see datastructures.
%        w    = angular frequency (default linspace(0,wc,257))
%        wc   = angular cutoff frequency (default 32/Tp)
%        data = [Hm0 Tp]
%               Hm0 = significant wave height (default 7 (m))
%               Tp  = peak period (default 11 (sec))
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
%  PMSPEC return the Pierson-Moskowitz spectral density:
%
%     S(w) = 5*Hm0^2/(wp*wn^5)*exp(-5/4*wn^-4)
% where
%     wp = 2*pi/Tp   and    wn = w/wp 
%
%  This is a suitable model for fully developed sea, i.e. a sea state
%  where the wind has been blowing long enough over a sufficiently open
%  stretch of water, so that the high-frequency waves have reached an
%  equilibrium. In the part of the spectrum where the frequency is
%  greater than the peak frequency (w>wp), the energy distribution is
%  proportional to w^-5.
%  The spectrum is identical with ITTC (International Towing Tank 
%  Conference), ISSC (International Ship and Offshore Structures Congress) 
%  and Bretschneider, wave spectrum given Hm0 and Tm01. It is also identical
%  with JONSWAP when the peakedness factor, gamma, is one.
%  For this spectrum, the following relations exist between the mean
%  period Tm01 = 2*pi*m0/m1, the peak period Tp and the mean
%  zero-upcrossing period Tz:
%
%           Tm01 = 1.086*Tz, Tp = 1.408*Tz and Tp=1.2965*Tm01
%
% Example: 
% S = pmspec(1.5,[6.5 10]); plotspec(S)
%
% See also  jonswap, torsethaugen, simpson

% References:
%

% Tested on: matlab 6.0, 5.3
% History:
% revised pab jan 2007
% -replaced code with call to ggamspec
% Revised pab April 2005
% -Changed parameterization to a generalized gamma formulation
% Revised pab 20.10.2001
% - initialized w to avoid Reference to uninitialized variable.
% - updated help header
% Revised by jr 03.04.2001
% - added wc to input
% - updated information
% By pab 01.12.99

monitor=0;
w = [];
if nargin<3||isempty(plotflag)
  plotflag=0;
end
if nargin<2||isempty(sdata)
  sdata=[7 11];
end %

if nargin<1||isempty(w1), 
   wc = 32/sdata(2);
elseif length(w1)==1,
   wc = w1; 
else
   w = w1(:);
end
nw = 257;
if isempty(w),
   w = linspace(0,wc,nw).'; 
end


n    = length(w);
S1   = createspec;
S1.S = zeros(n,1);
S1.w = w;
S1.norm=0; % The spectrum is not normalized


Hm0 = sdata(1);

Tp      = sdata(2);
S1.note = sprintf('Pierson-Moskowitz, Hm0 = %g, Tp = %g', Hm0,Tp);

if monitor
  disp(S1.note);
end


if Hm0>0
   N = 5;
   M = 4;
   wp = 2*pi/Tp;
   wn = w/wp;
   S1.S = (Hm0/4)^2/wp * ggamspec(wn,N,M);
end


if plotflag
  plotspec(S1,plotflag);
end
