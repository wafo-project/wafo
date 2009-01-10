function S1=ohspec(w1,sdata,plotflag)
%OHSPEC Calculates (and plots) a Ochi-Hubble spectral density.
% 
% CALL:  S = ohspec(w,data,plotflag); 
%        S = ohspec(wc,data,plotflag);
%
%        S    = a struct containing the spectral density, see datastructures.
%        w    = angular frequency (default linspace(0,33/Tp,257))
%        wc   = angular cutoff frequency (default 33/Tp)
%        data = [Hm0 Tp L]
%               Hm0 = significant wave height (default 7 (m))
%               Tp  = peak period (default 11 (sec))
%               L   = spectral shape parameter (default 3)
%    plotflag = 0, do not plot the spectrum (default).
%               1, plot the spectrum.
%
%  The OH spectrum parameterization used is 
%
%     S(w) = B^L*Hm0^2/(4*gamma(L)*w^(4*L+1))*exp(-B/w^4)
%  where
%       B = (L+1/4)*(2*pi/Tp)^4
%
% 
% Example: 
%  S = ohspec(0.95); 
%
% See also  ochihubble, jonswap, torsethaugen
  
% References:
% Ochi, M.K. and Hubble, E.N. (1976)
% 'On six-parameter wave spectra.'
% In Proc. 15th Conf. Coastal Engng., Vol.1, pp301-328

% Tested on: matlab 6.0, 5.3
% History:
% revised pab jan 2007
% -replaced code with call to ggamspec
% Revised pab Apr2005
% Fixed bug:
% Revised pab Dec2004
%  Fixed bug: w was previously not initiliazed.   
% Revised jr 03.04.2001
% - added wc to input
% - updated information
% Revised pab 24.11.2000
% - fixed bug: wrong sign L-0.25 is replaced with L+0.25
% revised pab 16.02.2000
% by pab 01.12.99

monitor=0;
w = [];
if nargin<3||isempty(plotflag)
  plotflag=0;
end
% Old call
%if nargin<1|isempty(w)
%  w=linspace(0,3,257).';
%end
sdata2=[7 11 3]; % default values
if nargin<2||isempty(sdata)
else
 ns1=length(sdata);
 k=find(~isnan(sdata(1:min(ns1,3))));
 if any(k)
   sdata2(k)=sdata(k);
 end
end %
if nargin<1||isempty(w1),
   wc = 33/sdata2(2);
elseif length(w1)==1,
   wc = w1; 
else
   w = w1 ;
end
nw = 257;
if isempty(w), 
   w = linspace(0,wc,nw).'; 
end
n       = length(w);

S1      = createspec;
S1.S    = zeros(n,1);
S1.w    = w(:);
S1.norm = 0; % The spectrum is not normalized


Hm0     = sdata2(1);
Tp      = sdata2(2);
L       = sdata2(3);
S1.note = ['Ochi-Hubble, Hm0 = ' num2str(Hm0)  ', Tp = ' num2str(Tp), ...
	', L = ' num2str(L)];
if monitor
  disp(['Hm0, Tp      = ' num2str([Hm0 Tp])])
end


if Hm0>0
  wp = 2*pi/Tp;
  % New call pab jan 2007
  wn = w/wp;
  N = 4*L+1;
  M = 4;
  
  S1.S = (Hm0/4)^2/wp * ggamspec(wn,N,M);
%   if 0,
%     k  = find(w>0);  % avoid division by zero
%     %Old call
%     %S1.S(k)=0.25*((L-0.25)*wp^4)^L/gamma(L)*Hm0^2./(w(k).^(4*L+1)).*exp(-(L-0.25)*(wp./w(k)).^4);
% 
%     % New call pab 24.11.2000
%     B = (L+0.25);
%     S1.S(k)=0.25*(B*wp^4)^L/gamma(L)*Hm0^2./(w(k).^(4*L+1)).*exp(-B*(wp./w(k)).^4);
%  end
end


if plotflag
  plotspec(S1,plotflag)
end
