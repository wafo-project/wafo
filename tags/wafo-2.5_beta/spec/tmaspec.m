function S1=tmaspec(w1,sdata,h,plotflag)
%TMASPEC Calculates a JONSWAP spectral density for finite water depth
%
% CALL:  S = tmaspec(w,data,h,plotflag); 
%        S = tmaspec(wc,data,h,plotflag);
%
%   S    = a struct containing the spectral density See datastructures
%   w    = angular frequency (default linspace(0,33/Tp,257))
%   wc   = angular cutoff frequency (default 33/Tp)
%   data = [Hm0 Tp gamma sa sb A], where
%          Hm0   = significant wave height (default 7 (m))
%          Tp    = peak period (default 11 (sec))
%          gamma = peakedness factor determines the concentraton
%                  of the spectrum on the peak frequency,  1 <= gamma <= 7. 
%                  (default depending on Hm0, Tp, see below)
%          sa,sb = spectral width parameters (default 0.07 0.09)
%          A     = alpha, normalization factor, (default -1)
%                  A<0 : A calculated by integration so that int S dw =Hm0^2/16
%                  A==0: A = 5.061*Hm0^2/Tp^4*(1-0.287*log(gamma))  
%                  A>0 : A = A
%          h     = waterdepth (default 42)
%  plotflag = 0, do not plot the spectrum (default).
%             1, plot the  spectrum.
%
%  For zero values, NaN's or parameters not specified in DATA the
%  default values are used. 
%
%  The evaluated spectrum is
%         S(w) = Sj(w)*phi(w,h)
%    where 
%         Sj  = jonswap spectrum
%         phi = modification due to water depth
%
%  The concept is base on a similarity law, and its validity is verified
%  through analysis of 3 data sets from: TEXEL, MARSEN projects (North
%  Sea) and ARSLOE project (Duck, North Carolina, USA)  
%
% Example:  
%      S = tmaspec([],[0 0 1])   % Bretschneider spectrum Hm0=7, Tp=11
%      S = tmaspec(1.5,[0 0 1])  % The same, cut at wc = 1.5
% 
% See also  jonswap, phi1, pmspec, torsethaugen

% References:
% Buows, E., Gunther, H., Rosenthal, W., and Vincent, C.L. (1985)
% 'Similarity of the wind wave spectrum in finite depth water: 1 spectral form.' 
%  J. Geophys. Res., Vol 90, No. C1, pp 975-986
%
% Hasselman et al. (1973)
% Measurements of Wind-Wave Growth and Swell Decay during the Joint
% North Sea Project (JONSWAP). 
% Ergansungsheft, Reihe A(8), Nr. 12, deutschen Hydrografischen
% Zeitschrift.

% Tested on: matlab 6.0, 5.2
% History:
% revised jr 03.04.2001
% - added wc to input
% - updated information
% revised pab 22.09.2000, corrected water depth, i.e. changed 
%    S1.h=inf to S1.h=h
% revised es 25.05.00 - some modifications to help text    
% by pab 16.02.2000

 
%NOTE: In order to calculate the short term statistics of the response,
%      it is extremely important that the resolution of the transfer
%      function is sufficiently good. In addition, the transfer function
%      must cover a sufficietn range of wave periods, especially in the
%      range where the wave spectrum contains most of its
%      energy. VIOLATION OF THIS MAY LEAD TO MEANINGLESS RESULTS FROM THE 
%      CALCULATIONS OF SHORT TERM STATISTICS. The highest wave period
%      should therefore be at least 2.5 to 3 times the highest peak
%      period in the transfer function. The lowest period should be selected 
%      so that the transfer function value is low. This low range is 
%      especially important when studying velocities and accelerations.

%monitor=0;

if nargin<4||isempty(plotflag)
  plotflag=0;
end
if nargin<3||isempty(h)
  h=42;
end



Hm0=7;Tp=11; gam=0; sa=0.07; sb=0.09; A=-1;% default values
data2=[Hm0 Tp gam sa sb A];
nd2=length(data2);
if (nargin>1) && ~isempty(sdata), 
  nd=length(sdata); 
  ind=find(~isnan(sdata(1:min(nd,nd2))));
  if any(ind) % replace default values with those from input data
    data2(ind)=sdata(ind);
  end
end
if (nd2>0) && (data2(1)>0), Hm0 = data2(1);end
if (nd2>1) && (data2(2)>0),  Tp = data2(2);end

w = [];
if nargin<1||isempty(w1), 
  wc = 33/Tp;
elseif length(w1)==1,
  wc = w1; 
else
  w = w1 ; 
end
nw = 257;
if isempty(w), w = linspace(0,wc,nw).'; end

% Old call
%if nargin<1|isempty(w)
%  w=linspace(0,33/Tp,257).';
%end

S1=jonswap(w,data2);
S1.h=h;
S1.norm=0; % The spectrum is not normalized
S1.note=['TMASPEC, Hm0 = ' num2str(Hm0)  ', Tp = ' num2str(Tp) ,', h = ' num2str(h)];
S1.S=S1.S.*phi1(S1.w,h);

return







function ph=phi1(w,h)
%PHI1 factor for transforming spectra to finite water depth spectra
%
% CALL: tr = phi1(w,h)
%
%     w = angular frequency
%     h = water depth
%
% Example: Transform a JONSWAP spectrum to a spectrum for waterdepth = 30 m
%   S = jonswap;
%   S1=S; S1.S=S1.S.*phi1(S1.w,30);

% reference
% Buows, E., Gunther, H., Rosenthal, W. and Vincent, C.L. (1985)
% 'Similarity of the wind wave spectrum in finite depth water: 1 spectral form.' 
%  J. Geophys. Res., Vol 90, No. C1, pp 975-986

% History:
% by pab 16.02.2000
g=gravity;
if h==inf,
  ph=ones(size(w));
  return
end
ph=zeros(size(w));

k1=w2k(w,0,inf);
dw1=2*w/g; % dw/dk|h=inf
k2=w2k(w,0,h);

dw2=ph;
ix=find(k1~=0);

dw2(ix)=dw1(ix)./(tanh(k2(ix)*h)+k2(ix)*h./cosh(k2(ix)*h).^2); % % dw/dk|h=h0
ph(ix)=(k1(ix)./k2(ix)).^3.*dw2(ix)./dw1(ix);



